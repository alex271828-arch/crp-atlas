"""Radius-only mode must never enter the metric/Gram search route."""
import builtins
import contextlib
import io
import json
from pathlib import Path
import shlex
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import sympy as sp

import benchmark
import run_body
import solve


class RadiusOnlyRouteTests(unittest.TestCase):
    @staticmethod
    def candidate():
        alpha = sp.Symbol('alpha', real=True)
        return {'status': 'ALGEBRAIC_CANDIDATE', 'derivation': [],
                'candidates': [{'primitive_symbol': alpha,
                                'primitive_polynomial': alpha,
                                'isolating_interval': (0, 0),
                                'squared_radii': [sp.Rational(3, 8)]}]}

    @staticmethod
    def unresolved():
        return {'status': 'UNRESOLVED', 'candidates': [],
                'reason': 'Mock radius search did not find an isolated point'}

    def exercise(self, cheap=None, radial=None, verification=(), *,
                 clearance=.1, equations_only=False):
        case = {'name': 'radius-only routing test', 'sources': {}, 'notes': [],
                'faces': [[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]],
                'groups': [[0, 1, 2, 3]], 'squared_radii_numeric': [.375],
                'preflight': {'center_clearance': clearance, 'face_pair_error': 0.}}
        model = {'radii_symbols': (sp.Symbol('s0', positive=True),), 'equations': []}
        reports = [{'verified': good, 'certified': good,
                    'status': 'CERTIFIED' if good else 'UNVERIFIED',
                    'reason': 'Mock independent verification'}
                   for good in verification]
        route_calls, geometry_imports = [], []
        real_import = builtins.__import__

        def disallow_geometry(name, *args, **kwargs):
            if name in ('gram_solve', 'gram_model'):
                geometry_imports.append(name)
                raise AssertionError('Radius-only mode must not import geometry search')
            return real_import(name, *args, **kwargs)

        def radial_search(model, **kwargs):
            is_cheap = kwargs.get('cheap_only', False)
            route_calls.append('cheap' if is_cheap else 'radial')
            return cheap if is_cheap else radial

        with tempfile.TemporaryDirectory() as folder:
            args = SimpleNamespace(out=folder, case='mock-case', bi=None,
                                   spheres=None, tolerance=1e-6, time_limit=120.,
                                   equations_only=equations_only, radius_only=True)
            with (patch('inputs.load_case', return_value=case),
                  patch('encoder.build_system', return_value=model) as encode,
                  patch('patches.enrich_system'),
                  patch('reduce.solve_model', side_effect=radial_search),
                  patch('verify.verify_solution', side_effect=reports) as verifier,
                  patch('builtins.__import__', side_effect=disallow_geometry),
                  contextlib.redirect_stdout(io.StringIO())):
                code = solve.worker(args)
                result = json.loads((Path(folder)/'result.json').read_text())
                files = {path.name for path in Path(folder).iterdir()}
                saved_reports = (json.loads((Path(folder)/'verification.json').read_text())
                                 if 'verification.json' in files else [])
                self.assertTrue(result['radius_only'])
                self.assertEqual(geometry_imports, [])
                self.assertNotIn('gram-derivation.json', files)
                return code, result, files, saved_reports, route_calls, encode.call_count, verifier.call_count

    def test_local_miss_proceeds_directly_to_full_radius_search(self):
        code, result, files, reports, calls, encodes, checks = self.exercise(
            self.unresolved(), self.candidate(), [True])
        self.assertEqual(code, 0, result)
        self.assertEqual(result['status'], 'SOLVED')
        self.assertEqual(result['search_route'], 'radial_only')
        self.assertEqual(calls, ['cheap', 'radial'])
        self.assertEqual((encodes, checks), (1, 1))
        self.assertEqual([r['search_route'] for r in reports], ['radial_only'])
        self.assertIn('certificate.json', files)

    def test_rejected_local_candidate_still_proceeds_to_full_radius_search(self):
        code, result, files, reports, calls, encodes, checks = self.exercise(
            self.candidate(), self.candidate(), [False, True])
        self.assertEqual(code, 0, result)
        self.assertEqual(result['search_route'], 'radial_only')
        self.assertEqual(calls, ['cheap', 'radial'])
        self.assertEqual(checks, 2)
        self.assertEqual([r['search_route'] for r in reports], ['local_radius', 'radial_only'])
        self.assertIn('certificate.json', files)

    def test_verified_local_candidate_finishes_without_full_search(self):
        code, result, files, reports, calls, encodes, checks = self.exercise(
            self.candidate(), self.unresolved(), [True])
        self.assertEqual(code, 0, result)
        self.assertEqual(result['search_route'], 'local_radius')
        self.assertEqual(calls, ['cheap'])
        self.assertEqual(checks, 1)
        self.assertIn('certificate.json', files)

    def test_failed_radius_search_stays_unresolved(self):
        code, result, files, reports, calls, encodes, checks = self.exercise(
            self.unresolved(), self.unresolved())
        self.assertEqual(code, 2)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertFalse(result['certified'])
        self.assertEqual(calls, ['cheap', 'radial'])
        self.assertEqual(checks, 0)
        self.assertNotIn('certificate.json', files)

    def test_noninterior_centers_are_rejected_even_for_equations_only(self):
        for clearance in (0., -.1):
            for equations_only in (False, True):
                with self.subTest(clearance=clearance, equations_only=equations_only):
                    code, result, files, reports, calls, encodes, checks = self.exercise(
                        clearance=clearance, equations_only=equations_only)
                    self.assertEqual(code, 2)
                    self.assertEqual(result['status'], 'UNRESOLVED')
                    self.assertIn('Radius-only mode requires a strictly interior sphere center', result['reason'])
                    self.assertEqual((calls, encodes, checks), ([], 0, 0))
                    self.assertNotIn('equations.json', files)
                    self.assertNotIn('certificate.json', files)

    def test_interior_equations_only_generates_radial_equations(self):
        code, result, files, reports, calls, encodes, checks = self.exercise(equations_only=True)
        self.assertEqual(code, 0, result)
        self.assertEqual(result['status'], 'EQUATIONS_ONLY')
        self.assertEqual((calls, encodes, checks), ([], 1, 0))
        self.assertIn('equations.json', files)
        self.assertNotIn('certificate.json', files)


class RadiusOnlyCommandTests(unittest.TestCase):
    def test_solve_parent_forwards_flag_to_worker(self):
        for radius_only in (False, True):
            with self.subTest(radius_only=radius_only), tempfile.TemporaryDirectory() as folder:
                argv = ['solve.py', '--case', 'mock-case', '--out', folder]
                if radius_only:
                    argv.append('--radius-only')

                def completed_worker(timeout):
                    (Path(folder)/'result.json').write_text(json.dumps({
                        'status': 'SOLVED', 'radius_only': radius_only, 'reason': 'Mock result'}))
                    return 0

                with (patch('sys.argv', argv), patch('solve.subprocess.Popen') as process,
                      contextlib.redirect_stdout(io.StringIO())):
                    process.return_value.wait.side_effect = completed_worker
                    self.assertEqual(solve.main(), 0)
                    command = process.call_args.args[0]
                    self.assertIn('--worker', command)
                    self.assertEqual('--radius-only' in command, radius_only)

    def test_run_body_forwards_flag(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            (directory/'input.json').write_text('{}')
            bi = directory/'numerical.json'
            bi.write_text('{}')
            partition = directory/'partition.json'
            partition.write_text('{}')
            data = {'input_file': str(bi), 'census_index': 153,
                    'display_id': 'Hmock', 'name': 'Mock body'}
            output = io.StringIO()
            with (patch('sys.argv', ['run_body.py', '153', '--radius-only', '--show-command']),
                  patch('run_body.select_body', return_value=(partition, data)),
                  contextlib.redirect_stdout(output)):
                self.assertEqual(run_body.main(), 0)
            command = shlex.split(output.getvalue().splitlines()[-1])
            self.assertIn('--radius-only', command)
            self.assertEqual(command[command.index('--spheres')+1], str(partition))

    def test_benchmark_forwards_flag_to_every_solver_run(self):
        with tempfile.TemporaryDirectory() as folder:
            with (patch('sys.argv', ['benchmark.py', '--out', folder, '--radius-only']),
                  patch('benchmark.subprocess.run', return_value=SimpleNamespace(returncode=2)) as run,
                  contextlib.redirect_stdout(io.StringIO())):
                self.assertEqual(benchmark.main(), 2)
            self.assertEqual(run.call_count, 10)
            for call in run.call_args_list:
                self.assertIn('--radius-only', call.args[0])


if __name__ == '__main__':
    unittest.main()
