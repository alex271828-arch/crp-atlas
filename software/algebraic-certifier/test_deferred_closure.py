"""Exact one-parameter angle closure after earlier radius eliminations."""
import unittest

import sympy as sp

from encoder import _relation_consequences, materialize_relation
from reduce import _time_limit, solve_model


def unit_corner(cosine):
    """A corner with an exactly unit cosine denominator."""
    return dict(N=cosine, M=sp.S.One, S=1-cosine**2,
                D_left=sp.S.One, D_right=sp.S.One)


class DeferredThreeAngleTests(unittest.TestCase):
    def setUp(self):
        self.budget = _time_limit(15)
        self.budget.__enter__()
        self.addCleanup(self.budget.__exit__, None, None, None)

    def test_signed_three_plus_one_plus_one_constant_identities(self):
        cases = [((3, 1, 1), 1, (sp.pi/4, sp.pi/6, sp.pi/12)),
                 ((3, -1, -1), 0, (sp.pi/5, sp.pi/5, 2*sp.pi/5)),
                 ((-3, 1, 1), 0, (sp.pi/5, sp.pi/5, 2*sp.pi/5)),
                 ((3, -1, 1), 1, (sp.pi/4, sp.pi/6, 5*sp.pi/12)),
                 ((-3, 1, -1), -1, (sp.pi/4, sp.pi/6, 5*sp.pi/12))]
        for counts, rhs, angles in cases:
            with self.subTest(counts=counts, rhs=rhs):
                self.assertEqual(sum(c*a for c, a in zip(counts, angles)), rhs*sp.pi)
                equations, method, _ = _relation_consequences(
                    counts, sp.Integer(rhs), [unit_corner(sp.cos(a)) for a in angles])
                self.assertEqual(len(equations), 1)
                self.assertIn('three-angle cosine Gram', method)
                self.assertEqual(sp.simplify(equations[0]), 0)
                self.assertFalse(equations[0].has(sp.Float))

    def test_signed_univariate_identities_are_exact_polynomial_identities(self):
        t = sp.Symbol('t', real=True)
        # If t=cos(A), the other cosines describe 2A, pi-2A, A, or pi-A.
        cases = [((3, -1, -1), 0, (t, 2*t**2-1, t)),
                 ((3, 1, 1), 2, (t, 1-2*t**2, -t)),
                 ((3, -1, 1), 1, (t, 2*t**2-1, -t)),
                 ((-3, 1, -1), -1, (t, 2*t**2-1, -t))]
        for counts, rhs, cosines in cases:
            with self.subTest(counts=counts, rhs=rhs):
                equations, _, _ = _relation_consequences(
                    counts, sp.Integer(rhs), [unit_corner(c) for c in cosines])
                self.assertEqual(len(equations), 1)
                self.assertEqual(sp.expand(equations[0]), 0)
                self.assertFalse(equations[0].has(sp.Float))

    def test_three_variable_five_angle_relation_stays_deferred(self):
        x, y, z = sp.symbols('x y z', real=True)
        equations, method, _ = _relation_consequences(
            (3, 1, 1), sp.Integer(2), [unit_corner(c) for c in (x, y, z)])
        self.assertEqual(equations, [])
        self.assertIn('without dense expansion', method)

    def test_paired_bivariate_five_angle_closure_is_exact(self):
        x, y = sp.symbols('x y', real=True)
        # 3A + (pi-2A) + (pi-A) = 2pi, with y=cos(pi-2A).
        equations, method, _ = _relation_consequences(
            (3, 1, 1), sp.Integer(2), [unit_corner(c) for c in (x, y, -x)])
        self.assertEqual(len(equations), 1)
        self.assertIn('without squaring', method)
        self.assertFalse(sp.Poly(equations[0], x, y).is_zero)
        self.assertEqual(sp.expand(equations[0].subs(y, 1-2*x**2)), 0)
        self.assertFalse(equations[0].has(sp.Float))

    def test_bivariate_five_angle_closure_with_unpaired_edges_stays_deferred(self):
        x, y = sp.symbols('x y', positive=True)
        corners = [unit_corner(x), unit_corner(y),
                   dict(N=x+y, M=x, S=x-(x+y)**2, D_left=x, D_right=sp.S.One)]
        equations, method, _ = _relation_consequences((3, 1, 1), sp.Integer(2), corners)
        self.assertEqual(equations, [])
        self.assertIn('without dense expansion', method)

    def test_bivariate_multiplicity_bound_is_retained(self):
        x, y = sp.symbols('x y', real=True)
        equations, _, _ = _relation_consequences(
            (5, 1, 1), sp.Integer(2), [unit_corner(c) for c in (x, y, -x)])
        self.assertEqual(equations, [])

    def test_body153_family_has_exact_quartic_closure_factor(self):
        d = sp.Symbol('d', real=True)
        r = sp.sqrt(5)
        # These are the first three fan corners of J48, after the cheap
        # metric relations give s0=a and s1=a+d. No numerical fit is used.
        a = (5-r)*d**2/8+(3+r)/2
        edge = 4*a-1
        adjacent = edge+2*d-d**2
        corners = []
        for numerator, left, right in ((2*a-1, edge, edge),
                                       (2*a+d-1, edge, adjacent),
                                       ((1-r)*a+d-1, edge, adjacent)):
            norm = left*right
            corners.append(dict(N=numerator, M=norm, S=norm-numerator**2,
                                D_left=left, D_right=right))
        equations, method, _ = materialize_relation(
            dict(counts=(3, 1, 1), rhs_pi=sp.Integer(2)), corners)
        self.assertEqual(len(equations), 1)
        self.assertIn('without squaring', method)
        numerator, denominator = sp.cancel(equations[0]).as_numer_denom()
        polynomial = sp.Poly(numerator, d, extension=r)
        quartic = sp.Poly(d**4+(3+7*r/5)*d**2-sp.Rational(9, 5)-4*r/5,
                         d, extension=r)
        self.assertFalse(polynomial.is_zero)
        self.assertTrue(polynomial.rem(quartic).is_zero)
        # Clearing the denominator does not hide a pole on that component.
        self.assertEqual(sp.gcd(sp.Poly(denominator, d, extension=r), quartic).degree(), 0)
        self.assertFalse(equations[0].has(sp.Float))

    def test_solver_retries_deferred_template_with_derived_radius_family(self):
        x, y, z = sp.symbols('x y z', real=True)
        radii = (x, y, z)
        # Equal corner angles satisfy 3A+A+A=2pi. The intended component
        # has cos(A)=cos(2pi/5), while the radii differ by exactly one.
        cosine = (sp.sqrt(5)-1)/4
        expected = [2+cosine, 3+cosine, 4+cosine]
        corners = [unit_corner(v-offset) for v, offset in zip(radii, (2, 3, 4))]
        relation = dict(counts=(3, 1, 1), rhs_pi=sp.Integer(2))
        self.assertEqual(materialize_relation(relation, corners)[0], [])
        model = dict(radius_variables=radii, variables=radii,
                     witness={v: float(value) for v, value in zip(radii, expected)},
                     consequence_equations=[y-x-1, z-x-2],
                     corner_types=corners, deferred_relations=[relation])
        result = solve_model(model, max_seconds=10)
        self.assertEqual(result['status'], 'ALGEBRAIC_CANDIDATE', result.get('reason'))
        values = result['candidate']['algebraic_values']
        for variable, value in zip(radii, expected):
            self.assertEqual(sp.to_number_field(values[variable]-value).coeffs(), [0])
        closure = next(p for p in result['phases']
                       if p['name'] == 'deferred_local_angle_closure')
        self.assertGreater(closure['solutions'], 0)
        self.assertFalse(any(p['name'] == 'compact_signed_radical_relations'
                             for p in result['phases']))


if __name__ == '__main__':
    unittest.main()
