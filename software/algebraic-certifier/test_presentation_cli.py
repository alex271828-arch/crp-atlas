"""Presentation budgets and reports must preserve completed certificates."""
import json
import contextlib
import io
import shlex
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import re

from present import run_presentation
import solve
import run_body
import benchmark
from check import verify_saved
from presentation_check import verify_presentation


class PresentationReportTests(unittest.TestCase):
    def source(self, folder):
        source = Path(folder)/'source'
        source.mkdir()
        saved = {'primitive_symbol': 'alpha', 'primitive_polynomial': 'alpha',
                 'isolating_interval': ['0', '0'], 'squared_radii': ['3/8'],
                 'groups': [[0, 1, 2, 3]], 'Q_coefficients': ['1', '-3/8']}
        (source/'radii.json').write_text(json.dumps(saved))
        (source/'result.json').write_text(json.dumps(
            {'status': 'SOLVED', 'certified': True, 'name': 'Tetrahedron <test>'}))
        return source, saved

    def test_search_timeout_keeps_exact_baseline_and_source_untouched(self):
        with tempfile.TemporaryDirectory() as folder:
            source, saved = self.source(folder)
            before = {p.name: p.read_bytes() for p in source.iterdir()}
            destination = Path(folder)/'display'
            report = run_presentation(source, destination, time_limit=.001)
            self.assertEqual(report['presentation_status'], 'VERIFIED')
            self.assertTrue(report['presentation_time_limit_reached'])
            self.assertEqual(before, {p.name: p.read_bytes() for p in source.iterdir()})
            display = json.loads((destination/'presentation.json').read_text())
            self.assertTrue(verify_presentation(saved, display)['verified'])
            page = (destination/'presentation.html').read_text()
            self.assertIn('<math ', page)
            self.assertIn('&lt;test&gt;', page)
            self.assertNotIn('<script', page)
            self.assertNotIn('<mfenced', page)
            for fragment in re.findall(r'<math\b.*?</math>', page, flags=re.S):
                ET.fromstring(fragment)
            self.assertIn('3', (destination/'presentation.md').read_text())

    def test_unresolved_source_cannot_be_presented_as_certified(self):
        with tempfile.TemporaryDirectory() as folder:
            source, _ = self.source(folder)
            (source/'result.json').write_text(json.dumps(
                {'status': 'UNRESOLVED', 'certified': False}))
            with self.assertRaisesRegex(ValueError, 'certified source'):
                run_presentation(source, Path(folder)/'display', time_limit=.001)
            self.assertFalse((Path(folder)/'display').exists())

    def test_unsafe_source_symbol_cannot_reach_mathml_on_timeout(self):
        with tempfile.TemporaryDirectory() as folder:
            source, saved = self.source(folder)
            saved['primitive_symbol'] = 'x</mi><script>alert(1)</script><mi>'
            (source/'radii.json').write_text(json.dumps(saved))
            destination = Path(folder)/'display'
            with self.assertRaisesRegex(ValueError, 'ASCII arithmetic identifier'):
                run_presentation(source, destination, time_limit=.001)
            self.assertFalse((destination/'presentation.html').exists())


class PresentationIntegrationTests(unittest.TestCase):
    @staticmethod
    def invoke(folder, flags=(), presentation=None):
        output = Path(folder)
        solved = {'status': 'SOLVED', 'certified': True, 'reason': 'Exact proof saved',
                  'seconds': .1, 'radius_only': True}

        def completed_worker(timeout):
            (output/'result.json').write_text(json.dumps(solved))
            (output/'certificate.json').write_text('Exact proof sentinel')
            return 0

        with (patch('sys.argv', ['solve.py', '--case', 'mock-case', '--out', folder,
                                 '--radius-only', *flags]),
              patch('solve.subprocess.Popen') as process,
              patch('present.run_presentation', side_effect=presentation) as display,
              contextlib.redirect_stdout(io.StringIO())):
            process.return_value.wait.side_effect = completed_worker
            code = solve.main()
        return code, json.loads((output/'result.json').read_text()), display, process

    def test_optional_presentation_failure_preserves_solved_certificate(self):
        with tempfile.TemporaryDirectory() as folder:
            code, result, display, process = self.invoke(
                folder, presentation=RuntimeError('Search unavailable'))
            self.assertEqual(code, 0)
            self.assertEqual(result['status'], 'SOLVED')
            self.assertTrue(result['certified'])
            self.assertEqual(result['seconds'], .1)
            self.assertEqual(result['presentation_status'], 'UNAVAILABLE')
            self.assertEqual((Path(folder)/'certificate.json').read_text(), 'Exact proof sentinel')
            display.assert_called_once_with(Path(folder).resolve(), time_limit=5.)
            process.return_value.wait.assert_called_once_with(timeout=120.)

    def test_presentation_has_its_own_budget_and_keeps_timeout_result_solved(self):
        with tempfile.TemporaryDirectory() as folder:
            report = {'presentation_status': 'VERIFIED', 'presentation_improved': False,
                      'presentation_time_limit_reached': True, 'presentation_seconds': 7.}
            code, result, display, process = self.invoke(
                folder, ['--presentation-time-limit', '7'], presentation=lambda *a, **kw: report)
            self.assertEqual(code, 0)
            self.assertEqual(result['status'], 'SOLVED')
            self.assertTrue(result['certified'])
            self.assertTrue(result['presentation_time_limit_reached'])
            display.assert_called_once_with(Path(folder).resolve(), time_limit=7.)
            process.return_value.wait.assert_called_once_with(timeout=120.)

    def test_interrupted_presentation_keeps_completed_proof_solved(self):
        with tempfile.TemporaryDirectory() as folder:
            code, result, display, process = self.invoke(folder, presentation=KeyboardInterrupt())
            self.assertEqual(code, 130)
            self.assertEqual(result['status'], 'SOLVED')
            self.assertTrue(result['certified'])
            self.assertEqual(result['presentation_status'], 'INTERRUPTED')

    def test_disabled_presentation_and_overwrite_remove_old_display_artifacts(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            (output/'result.json').write_text('{}')
            stale = ('presentation.json', 'presentation.md', 'presentation.html',
                     'presentation.log', 'presentation.json.tmp')
            for name in stale:
                (output/name).write_text('old display')
            code, result, display, process = self.invoke(folder, ['--overwrite', '--no-presentation'])
            self.assertEqual(code, 0)
            self.assertEqual(result['status'], 'SOLVED')
            display.assert_not_called()
            self.assertFalse(any((output/name).exists() for name in stale))

    def test_run_body_forwards_both_presentation_controls(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            (directory/'input.json').write_text('{}')
            bi = directory/'numerical.json'
            bi.write_text('{}')
            partition = directory/'partition.json'
            partition.write_text('{}')
            data = {'input_file': str(bi), 'census_index': 1,
                    'display_id': 'Hmock', 'name': 'Mock body'}
            output = io.StringIO()
            with (patch('sys.argv', ['run_body.py', '1', '--radius-only', '--no-presentation',
                                     '--presentation-time-limit', '9', '--show-command']),
                  patch('run_body.select_body', return_value=(partition, data)),
                  contextlib.redirect_stdout(output)):
                self.assertEqual(run_body.main(), 0)
            command = shlex.split(output.getvalue().splitlines()[-1])
            self.assertIn('--no-presentation', command)
            self.assertIn('--radius-only', command)
            self.assertEqual(command[command.index('--presentation-time-limit')+1], '9.0')

    def test_benchmark_forwards_both_presentation_controls(self):
        with tempfile.TemporaryDirectory() as folder:
            with (patch('sys.argv', ['benchmark.py', '--out', folder, '--radius-only',
                                     '--no-presentation', '--presentation-time-limit', '9']),
                  patch('benchmark.subprocess.run', return_value=SimpleNamespace(returncode=2)) as run,
                  contextlib.redirect_stdout(io.StringIO())):
                self.assertEqual(benchmark.main(), 2)
            self.assertEqual(run.call_count, 10)
            for call in run.call_args_list:
                command = call.args[0]
                self.assertIn('--no-presentation', command)
                self.assertIn('--radius-only', command)
                self.assertEqual(command[command.index('--presentation-time-limit')+1], '9.0')

    def test_saved_checker_accepts_identity_and_rejects_tampered_display(self):
        from present import baseline
        with tempfile.TemporaryDirectory() as folder:
            source, saved = PresentationReportTests().source(folder)
            (source/'input.json').write_text(json.dumps({
                'faces': [[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]],
                'groups': [[0, 1, 2, 3]]}))
            self.assertTrue(verify_saved(source)['verified'])
            display = baseline(saved)
            (source/'presentation.json').write_text(json.dumps(display))
            report = verify_saved(source)
            self.assertTrue(report['verified'], report)
            self.assertTrue(report['presentation_check']['verified'])
            display['squared_radii'][0] = '1/2'
            (source/'presentation.json').write_text(json.dumps(display))
            report = verify_saved(source)
            self.assertFalse(report['verified'], report)
            self.assertEqual(report['status'], 'UNVERIFIED')
            self.assertFalse(report['presentation_check']['verified'])


if __name__ == '__main__':
    unittest.main()
