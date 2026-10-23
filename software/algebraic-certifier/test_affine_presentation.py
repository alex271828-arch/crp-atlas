"""Readable affine radius families must retain exact field identities."""
import unittest

import sympy as sp

from presentation import (_anchored_forms, _check_candidate, _field_value,
                          _score, build_presentation)
from presentation_check import verify_presentation


class AffinePresentationTests(unittest.TestCase):
    def setUp(self):
        self.t = sp.Symbol('t')
        t = self.t
        self.poly = sp.Poly(25*t**8+150*t**6-110*t**4+10*t**2+1, t, domain=sp.QQ)
        self.a = (5*t**6+30*t**4-5*t**2-12)/(8*(7*t**2-4))
        self.b = (-5*t**4-8*t**2+5)/(14*t**2-8)

    def radius_values(self, include_fourth=False):
        t, a, b = self.t, self.a, self.b
        values = [a, a+t, a+b*t]
        if include_fourth:
            values.append(a+(2*b-1)*t+1)
        return [_field_value(value, t, self.poly) for value in values]

    def test_anchor_recovers_shared_offset_and_generator_difference(self):
        values = self.radius_values()
        definitions, forms = _anchored_forms(self.poly, values)
        self.assertEqual(forms[0], sp.Symbol('a'))
        self.assertEqual(sp.expand(forms[1]-forms[0]), self.t)
        self.assertLess(_score(self.poly.as_expr(), definitions, forms), 160)
        _check_candidate(self.poly, values, sp.Poly(self.t, self.t),
                         self.poly, definitions, forms)

    def test_fourth_radius_shares_an_affine_quotient_and_constant(self):
        values = self.radius_values(include_fourth=True)
        definitions, forms = _anchored_forms(self.poly, values)
        self.assertEqual(sp.expand(forms[3]-2*forms[2]+forms[1]-1), 0)
        self.assertLess(_score(self.poly.as_expr(), definitions, forms), 200)
        _check_candidate(self.poly, values, sp.Poly(self.t, self.t),
                         self.poly, definitions, forms)

    def test_every_anchor_and_rational_field_remain_exact(self):
        values = self.radius_values(include_fourth=True)
        for anchor in range(len(values)):
            with self.subTest(anchor=anchor):
                definitions, forms = _anchored_forms(self.poly, values, anchor)
                _check_candidate(self.poly, values, sp.Poly(self.t, self.t),
                                 self.poly, definitions, forms)
        poly = sp.Poly(self.t, self.t, domain=sp.QQ)
        values = [sp.Poly(v, self.t, domain=sp.QQ) for v in
                  [sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(11, 6)]]
        definitions, forms = _anchored_forms(poly, values)
        _check_candidate(poly, values, sp.Poly(0, self.t), poly, definitions, forms)

    def test_complete_engine_and_independent_checker_accept_short_family(self):
        alpha = sp.Symbol('alpha')
        saved = {'primitive_symbol': 'alpha',
                 'primitive_polynomial': str(self.poly.as_expr().subs(self.t, alpha)),
                 'isolating_interval': ['733/1000', '734/1000'],
                 'squared_radii': [str(v.as_expr().subs(self.t, alpha))
                                   for v in self.radius_values()]}
        checkpoints = []
        result = build_presentation(saved, max_seconds=3, on_improvement=checkpoints.append)
        self.assertTrue(verify_presentation(saved, result)['verified'])
        self.assertLess(result['score'], 160)
        self.assertTrue(all(a['score'] > b['score'] for a, b in
                            zip(checkpoints, checkpoints[1:])))


if __name__ == '__main__':
    unittest.main()
