"""Two-radius fan closure, using only radial equations and corner cosines."""
import unittest

import sympy as sp

from encoder import materialize_relation, _relation_consequences
from reduce import solve_model


def corner(numerator, left, right):
    norm = left*right
    return dict(N=numerator, M=norm, S=norm-numerator**2,
                D_left=left, D_right=right)


def difference_polynomial(d):
    return 32*d**8 + 128*d**5 - 48*d**4 - 224*d**3 - 80*d**2 + 1


def radial_model(*, extra_radius=False):
    """The two independent fan equations of the square bicupola example.

    Numerical radii select a branch only.  All equation coefficients follow
    from triangle and square corner distances.  The optional third radius
    introduces an exact affine substitution before the same fan is usable.
    """
    radii = sp.symbols('s0:3' if extra_radius else 's0:2', positive=True)
    a, b = radii[:2]
    d = b-a
    first_fan = (16*a**2 - 16*a*d**2 + 32*a*d - 40*a
                 + 4*d**4 - 16*d**3 + 28*d**2 - 24*d + 17)
    across = 4*a*b-(a+b-1)**2
    outside = 4*b-1
    third = radii[2]-1 if extra_radius else b
    corners = [corner(a+b-1, across, outside),
               corner(a-b-1, across, outside),
               corner(2*third-1, 4*third-1, 4*third-1)]
    equations = [sp.expand(first_fan)]
    numeric = [1.7933479490264195, 1.8921338988540684]
    if extra_radius:
        equations.append(radii[2]-b-1)
        numeric.append(numeric[1]+1)
    return dict(radius_variables=radii, variables=radii,
                witness=dict(zip(radii, numeric)),
                consequence_equations=equations, corner_types=corners,
                deferred_relations=[dict(counts=(1, 1, 3), rhs_pi=sp.Integer(2))])


class BivariateRadialClosureTests(unittest.TestCase):
    def assert_radial_solution(self, model):
        result = solve_model(model, max_seconds=10)
        self.assertEqual(result['status'], 'ALGEBRAIC_CANDIDATE', result.get('reason'))
        candidate = result['candidate']
        alpha = candidate['primitive_symbol']
        polynomial = sp.Poly(candidate['primitive_polynomial'], alpha, domain=sp.QQ)
        self.assertEqual(polynomial.monic(),
                         sp.Poly(difference_polynomial(alpha), alpha).monic())
        values = candidate['radius_expressions']
        a, b = model['radius_variables'][:2]
        self.assertTrue(sp.Poly(values[b]-values[a]-alpha, alpha).rem(polynomial).is_zero)
        closure = next(phase for phase in result['phases']
                       if phase['name'] == 'deferred_local_angle_closure')
        self.assertGreater(closure['solutions'], 0)
        self.assertFalse(any(phase['name'] == 'compact_signed_radical_relations'
                             for phase in result['phases']))
        equations = list(model['consequence_equations'])
        for relation in model['deferred_relations']:
            # Recover the original fan after the exact affine substitution,
            # if one was required to bring it inside the construction bound.
            substitutions = {}
            if len(model['radius_variables']) == 3:
                substitutions[model['radius_variables'][2]] = b+1
            equations.extend(materialize_relation(
                relation, model['corner_types'], substitutions)[0])
        for equation in equations:
            numerator, denominator = sp.cancel(equation.subs(values)).as_numer_denom()
            self.assertTrue(sp.Poly(numerator, alpha, domain=sp.QQ).rem(polynomial).is_zero)
            self.assertEqual(sp.gcd(sp.Poly(denominator, alpha, domain=sp.QQ),
                                    polynomial).degree(), 0)
        return result

    def test_two_radius_fans_isolate_the_exact_difference_field(self):
        self.assert_radial_solution(radial_model())

    def test_affine_reduction_unlocks_a_deferred_bivariate_fan(self):
        model = radial_model(extra_radius=True)
        relation = model['deferred_relations'][0]
        self.assertEqual(materialize_relation(relation, model['corner_types'])[0], [])
        result = self.assert_radial_solution(model)
        self.assertTrue(any(entry.get('substituted_radius_family_first') is True
                            for entry in result['derivation']
                            if entry['operation'] == 'local_angle_template_instantiation'))

    def test_signed_identities_with_positive_variable_denominator_scales(self):
        x = sp.Symbol('x', real=True)
        scale = sp.Symbol('scale', positive=True)
        scales = (scale, scale+1, scale**2+1)
        cases = [((3, -1, -1), 0, (x, 2*x**2-1, x)),
                 ((3, -1, 1), 1, (x, 2*x**2-1, -x)),
                 ((-3, 1, -1), -1, (x, 2*x**2-1, -x))]
        for counts, rhs, cosines in cases:
            with self.subTest(counts=counts, rhs=rhs):
                corners = [corner(c*q, q, q) for c, q in zip(cosines, scales)]
                equations, method, _ = _relation_consequences(counts, sp.Integer(rhs), corners)
                self.assertEqual(len(equations), 1)
                self.assertIn('without squaring', method)
                self.assertEqual(sp.cancel(equations[0]), 0)
                self.assertFalse(equations[0].has(sp.Float))


if __name__ == '__main__':
    unittest.main()
