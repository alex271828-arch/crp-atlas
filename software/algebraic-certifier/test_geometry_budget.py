"""Geometry uses the available search budget when radial fallback is unavailable."""
import contextlib
import io
import itertools
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import solve


class GeometryBudgetTests(unittest.TestCase):
    def exercise_worker(self, clearance, time_limit, elapsed=0.):
        case = {'name': 'geometry budget test', 'sources': {}, 'notes': [],
                'faces': [], 'groups': [], 'squared_radii_numeric': [],
                'preflight': {'center_clearance': clearance, 'face_pair_error': 0.}}
        unresolved = {'status': 'UNRESOLVED', 'candidates': [],
                      'reason': 'No mock candidate'}
        progress_records = []
        original_progress = solve.progress

        def record_progress(out, stage, **details):
            progress_records.append(dict(stage=stage, **details))
            original_progress(out, stage, **details)

        with tempfile.TemporaryDirectory() as folder:
            args = SimpleNamespace(out=folder, case='mock-case', bi=None,
                                   spheres=None, tolerance=1e-6, time_limit=time_limit,
                                   equations_only=False, radius_only=False)
            with (patch('inputs.load_case', return_value=case),
                  patch('encoder.build_system', return_value={}),
                  patch('patches.enrich_system'),
                  patch('reduce._time_limit', side_effect=lambda _: contextlib.nullcontext()),
                  patch('reduce.solve_model', return_value=unresolved) as radial,
                  patch('gram_solve.solve_gram_model', return_value=unresolved) as geometry,
                  patch('solve.certify_candidates', return_value=False),
                  patch('solve.time.monotonic', side_effect=itertools.chain(
                      [100.], itertools.repeat(100.+elapsed))),
                  patch('solve.progress', side_effect=record_progress),
                  contextlib.redirect_stdout(io.StringIO())):
                self.assertEqual(solve.worker(args), 2)
                result = json.loads((Path(folder)/'result.json').read_text())
                self.assertEqual(result['status'], 'UNRESOLVED')
                self.assertFalse(result['certified'])
                geometry.assert_called_once()
                budget = geometry.call_args.kwargs['max_seconds']
                stage = next(p for p in progress_records
                             if p['stage'] == 'solving_exact_distance_and_symmetry_equations')
                self.assertEqual(stage['search_budget_seconds'], budget)
                return budget, stage['radial_fallback_available'], radial.call_args_list

    def test_noninterior_centers_use_available_time_without_sixty_five_second_cap(self):
        for clearance in (0., -.1, 1e-8):
            for time_limit in (20., 120., 600.):
                with self.subTest(clearance=clearance, time_limit=time_limit):
                    budget, radial_available, radial_calls = self.exercise_worker(clearance, time_limit)
                    self.assertEqual(budget, time_limit-2.)
                    self.assertFalse(radial_available)
                    self.assertEqual(radial_calls, [])

    def test_geometry_budget_excludes_time_already_spent(self):
        budget, radial_available, radial_calls = self.exercise_worker(0., 120., elapsed=17.)
        self.assertEqual(budget, 101.)
        self.assertFalse(radial_available)
        self.assertEqual(radial_calls, [])

    def test_interior_centers_retain_cap_and_fraction_for_radial_fallback(self):
        for time_limit, expected in ((20., 11.7), (120., 65.), (600., 65.)):
            with self.subTest(time_limit=time_limit):
                budget, radial_available, radial_calls = self.exercise_worker(.1, time_limit)
                self.assertAlmostEqual(budget, expected)
                self.assertTrue(radial_available)
                self.assertEqual(len(radial_calls), 2)
                self.assertTrue(radial_calls[0].kwargs['cheap_only'])
                self.assertNotIn('cheap_only', radial_calls[1].kwargs)

    def test_parent_still_forwards_and_enforces_requested_hard_limit(self):
        for time_limit in (20., 120., 600.):
            with self.subTest(time_limit=time_limit), tempfile.TemporaryDirectory() as folder:
                argv = ['solve.py', '--case', 'mock-case', '--out', folder,
                        '--time-limit', str(time_limit)]

                def completed_worker(timeout):
                    (Path(folder)/'result.json').write_text(json.dumps({
                        'status': 'UNRESOLVED', 'certified': False, 'reason': 'Mock result'}))
                    return 2

                with (patch('sys.argv', argv), patch('solve.subprocess.Popen') as process,
                      contextlib.redirect_stdout(io.StringIO())):
                    process.return_value.wait.side_effect = completed_worker
                    self.assertEqual(solve.main(), 2)
                    command = process.call_args.args[0]
                    self.assertEqual(command[command.index('--time-limit')+1], str(time_limit))
                    process.return_value.wait.assert_called_once_with(timeout=time_limit)


if __name__ == '__main__':
    unittest.main()
