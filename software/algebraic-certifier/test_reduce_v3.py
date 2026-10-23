"""Generic local metric reductions and lazy angle encoding."""
import unittest
from pathlib import Path
from unittest.mock import patch
import sympy as sp
from encoder import build_system, materialize_relation
from inputs import load_case
from reduce import _staged_linear, _substitute, _select_factors, solve_low_degree


class LocalMetricReductionTests(unittest.TestCase):
    def test_two_sphere_intersections_produce_exact_coordinates(self):
        x, y = sp.symbols('x y', real=True)
        result = solve_low_degree([x*x+y*y-1, (x-1)**2+y*y-1], [x, y],
                                  {x: .5, y: 3**.5/2}, max_seconds=3)
        self.assertEqual(result['status'], 'ALGEBRAIC_CANDIDATE', result.get('reason'))
        solution = result['solutions'][0]
        self.assertEqual(solution[x], sp.Rational(1, 2))
        self.assertEqual(sp.to_number_field(solution[y]**2).as_expr(), sp.Rational(3, 4))
        self.assertTrue(any(row['operation'] == 'quadratic_leading_term_cancellation'
                            for row in result['derivation']))

    def test_local_quadratics_close_without_a_resultant(self):
        u, z = sp.symbols('u z', real=True)
        h = 1-u*u
        equations = [h*(z*z-2*z)-u+sp.Rational(1, 2),
                     4*h*((2*u*u-1)*z-u)**2+4*u*u-3]
        guide = {u: .8527269428464168}
        guide[z] = 1+(1+(guide[u]-.5)/(1-guide[u]**2))**.5
        with patch.object(sp, 'resultant', side_effect=AssertionError('no resultants')):
            result = solve_low_degree(equations, [u, z], guide, max_seconds=5)
        self.assertEqual(result['status'], 'ALGEBRAIC_CANDIDATE', result.get('reason'))
        roots = result['solutions'][0][u].atoms(sp.CRootOf)
        self.assertEqual(len(roots), 1)
        self.assertEqual(next(iter(roots)).poly.degree(), 4)
        self.assertTrue(any(row['operation'] == 'linear_elimination' and row['variable'] == z
                            for row in result['derivation']))

    def test_partial_family_is_preserved(self):
        x, y, z = sp.symbols('x y z')
        result = solve_low_degree([x-y-z], [x, y, z], {x: 3., y: 1., z: 2.},
                                  max_seconds=2)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertEqual(result['solutions'], [])
        self.assertEqual(len(result['remaining_variables']), 2)
        self.assertEqual(len(result['substitutions']), 1)


class DeferredRadiusRelationsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        case_path = Path(__file__).parent/'examples'/'J88'
        case = load_case(case_path, spheres=case_path/'spheres.json')
        cls.model = build_system(case['faces'], case['groups'],
                                 case['squared_radii_numeric'], defer_dense=True)

    def test_short_relations_derive_radius_family_before_dense_closure(self):
        model = self.model
        radii = model['radii_symbols']
        differences = sp.symbols('dr1:4', real=True)
        transform = {radii[i]: radii[0]+sum(differences[:i]) for i in range(1, 4)}
        guide = dict(model['witness'])
        guide.update({differences[i-1]: guide[radii[i]]-guide[radii[i-1]]
                      for i in range(1, 4)})
        trace = []
        equations = [_substitute(_select_factors(e, guide, trace), transform)
                     for e in model['consequence_equations']]
        equations, remaining, substitutions = _staged_linear(
            equations, [radii[0], *differences], guide, trace)
        family = {v: sp.cancel(transform.get(v, v).xreplace(substitutions)) for v in radii}
        e = differences[-1]
        expected_a = (e**4+3*e**2-1)/(4*e**2-1)
        expected_c = expected_a+(1-2*e**2)*(1+e)/(4*e**2-1)
        self.assertEqual(sp.cancel(family[radii[0]]-expected_a), 0)
        self.assertEqual(sp.cancel(family[radii[2]]-expected_c), 0)
        self.assertEqual(sp.cancel(family[radii[3]]-expected_c-e), 0)
        self.assertEqual(equations, [])
        self.assertEqual(len(remaining), 2)
        self.assertTrue(model['deferred_relations'])
        self.assertLessEqual(max(sp.total_degree(eq) for eq in model['consequence_equations']), 7)

    def test_deferred_materialization_preserves_original_relation(self):
        relation = next(r for r in self.model['deferred_relations']
                        if sorted(abs(c) for c in r['counts'] if c) == [1, 2, 2])
        substitution = {self.model['radii_symbols'][0]: sp.Rational(3, 5)}
        original, _, _ = materialize_relation(relation, self.model['corner_types'])
        delayed, _, _ = materialize_relation(relation, self.model['corner_types'], substitution)
        self.assertEqual(len(original), len(delayed))
        for a, b in zip(original, delayed):
            self.assertEqual(sp.cancel(a.xreplace(substitution)-b), 0)
            self.assertFalse(b.has(sp.Float))


if __name__ == '__main__':
    unittest.main()
