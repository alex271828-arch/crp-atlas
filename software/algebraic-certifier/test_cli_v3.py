"""Route acceptance requires independent verification, not algebraic status."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import sympy as sp

from solve import worker


class VerifiedRouteTests(unittest.TestCase):
    @staticmethod
    def candidate():
        alpha = sp.Symbol('alpha', real=True)
        return {'status': 'ALGEBRAIC_CANDIDATE',
                'candidates': [{'primitive_symbol': alpha,
                                'primitive_polynomial': alpha,
                                'isolating_interval': (0, 0),
                                'squared_radii': [sp.Rational(3, 8)]}],
                'derivation': []}

    @staticmethod
    def unresolved():
        return {'status': 'UNRESOLVED', 'candidates': [],
                'reason': 'Mock search did not find an isolated point'}

    def exercise(self, cheap, gram, radial, verification):
        case = {'name': 'routing test', 'sources': {}, 'notes': [],
                'faces': [[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]],
                'groups': [[0, 1, 2, 3]], 'squared_radii_numeric': [.375],
                'preflight': {'center_clearance': .1, 'face_pair_error': 0.}}
        model = {'radii_symbols': (sp.Symbol('s0', positive=True),), 'equations': []}
        reports = [{'verified': good, 'certified': good,
                    'status': 'CERTIFIED' if good else 'UNVERIFIED',
                    'reason': 'Mock independent verification'}
                   for good in verification]
        route_calls = []

        def radial_search(model, **kwargs):
            route_calls.append('cheap' if kwargs.get('cheap_only') else 'radial')
            return cheap if kwargs.get('cheap_only') else radial

        with tempfile.TemporaryDirectory() as folder:
            args = SimpleNamespace(out=folder, case='mock-case', bi=None,
                                   spheres=None, tolerance=1e-6,
                                   time_limit=120., equations_only=False)
            with (patch('inputs.load_case', return_value=case),
                  patch('encoder.build_system', return_value=model),
                  patch('patches.enrich_system'),
                  patch('reduce.solve_model', side_effect=radial_search),
                  patch('gram_solve.solve_gram_model', return_value=gram) as gram_search,
                  patch('verify.verify_solution', side_effect=reports) as verifier,
                  contextlib.redirect_stdout(io.StringIO())):
                code = worker(args)
                result = json.loads((Path(folder)/'result.json').read_text())
                saved_reports = (json.loads((Path(folder)/'verification.json').read_text())
                                 if (Path(folder)/'verification.json').exists() else [])
                certificate_exists = (Path(folder)/'certificate.json').exists()
                return code, result, saved_reports, certificate_exists, route_calls, gram_search.call_count, verifier.call_count

    def test_rejected_cheap_candidate_proceeds_to_gram(self):
        code, result, reports, certificate, calls, grams, checks = self.exercise(
            self.candidate(), self.candidate(), self.unresolved(), [False, True])
        self.assertEqual(code, 0)
        self.assertEqual(result['status'], 'SOLVED')
        self.assertEqual(result['search_route'], 'symmetry_gram')
        self.assertEqual(calls, ['cheap'])
        self.assertEqual((grams, checks), (1, 2))
        self.assertEqual([r['search_route'] for r in reports], ['local_radius', 'symmetry_gram'])
        self.assertTrue(certificate)

    def test_rejected_gram_candidate_proceeds_to_radial_fallback(self):
        code, result, reports, certificate, calls, grams, checks = self.exercise(
            self.unresolved(), self.candidate(), self.candidate(), [False, True])
        self.assertEqual(code, 0)
        self.assertEqual(result['status'], 'SOLVED')
        self.assertEqual(result['search_route'], 'radial_fallback')
        self.assertEqual(calls, ['cheap', 'radial'])
        self.assertEqual((grams, checks), (1, 2))
        self.assertEqual([r['search_route'] for r in reports], ['symmetry_gram', 'radial_fallback'])
        self.assertTrue(certificate)

    def test_all_rejected_candidates_remain_unresolved(self):
        code, result, reports, certificate, calls, grams, checks = self.exercise(
            self.candidate(), self.candidate(), self.candidate(), [False, False, False])
        self.assertEqual(code, 2)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertFalse(result['certified'])
        self.assertEqual(result['candidates_checked'], 3)
        self.assertEqual(calls, ['cheap', 'radial'])
        self.assertEqual((grams, checks), (1, 3))
        self.assertEqual(len(reports), 3)
        self.assertFalse(certificate)


class NoninteriorRouteTests(unittest.TestCase):
    def run_boundary_case(self, *, correct_radius):
        # A unit tetrahedron with its common center on the face 0,1,2.
        # Only the search is mocked; exact geometric certification is real.
        alpha = sp.Symbol('alpha', real=True)
        groups = [[0, 1, 2], [3]]
        case = {'name': 'boundary center routing test', 'sources': {}, 'notes': [],
                'faces': [[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]],
                'groups': groups, 'squared_radii_numeric': [1/3, 2/3],
                'preflight': {'center_clearance': 0., 'face_pair_error': 0.}}
        half = sp.Rational(1, 2)
        third = sp.Rational(1, 3)
        geometry = {'metric': [[1, half, half], [half, 1, half], [half, half, 1]],
                    'coordinates': [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    'center': [third, third, 0]}
        solution = {'status': 'ALGEBRAIC_CANDIDATE', 'derivation': [],
                    'candidates': [{'primitive_symbol': alpha, 'primitive_polynomial': alpha,
                                    'isolating_interval': (0, 0),
                                    'squared_radii': [third, 2*third if correct_radius else sp.Rational(3, 4)],
                                    'geometric_certificate': geometry}]}
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(out=directory, case='mock-case', bi=None, spheres=None,
                                   tolerance=1e-6, time_limit=120., equations_only=False)
            with (patch('inputs.load_case', return_value=case),
                  patch('encoder.build_system', side_effect=AssertionError('must bypass radial encoding')) as encode,
                  patch('reduce.solve_model', side_effect=AssertionError('must bypass radial solving')) as radial,
                  patch('verify.verify_solution', side_effect=AssertionError('must bypass radial certification')) as radial_check,
                  patch('gram_solve.solve_gram_model', return_value=solution) as gram,
                  contextlib.redirect_stdout(io.StringIO())):
                status = worker(args)
                result = json.loads((Path(directory)/'result.json').read_text())
                reports = json.loads((Path(directory)/'verification.json').read_text())
                saved = (json.loads((Path(directory)/'radii.json').read_text())
                         if (Path(directory)/'radii.json').exists() else None)
                certificate = (Path(directory)/'certificate.json').exists()
                encode.assert_not_called()
                radial.assert_not_called()
                radial_check.assert_not_called()
                gram.assert_called_once()
                return status, result, reports, saved, certificate

    def test_boundary_center_bypasses_radial_and_certifies_exact_geometry(self):
        status, result, reports, saved, certificate = self.run_boundary_case(correct_radius=True)
        self.assertEqual(status, 0, result)
        self.assertEqual(result['status'], 'SOLVED')
        self.assertTrue(result['certified'])
        self.assertEqual(result['verification_method'], 'exact_gram_geometry')
        self.assertEqual(reports[-1]['center_location'], 'boundary')
        self.assertEqual(saved['verification_method'], 'exact_gram_geometry')
        self.assertIn('geometric_certificate', saved)
        self.assertTrue(certificate)

    def test_invalid_boundary_center_candidate_remains_unresolved(self):
        status, result, reports, saved, certificate = self.run_boundary_case(correct_radius=False)
        self.assertEqual(status, 2)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertFalse(result['certified'])
        self.assertFalse(reports[-1]['verified'])
        self.assertIsNone(saved)
        self.assertFalse(certificate)


if __name__ == '__main__':
    unittest.main()
