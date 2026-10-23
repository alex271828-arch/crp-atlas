"""Fast algebraic controls for the generic radius reduction module.

Run with ``python -m unittest discover -s research/algebraic_radii_v3``.
The four real-body integrations are exercised separately by the benchmark.
"""
import unittest
from unittest.mock import patch

import sympy as sp

from reduce import _substitute, solve_model


class ReductionTests(unittest.TestCase):
    def setUp(self):
        self.x, self.y = sp.symbols("s0 s1", positive=True)

    def model(self, equations, witness):
        return dict(radius_variables=(self.x, self.y),
                    variables=(self.x, self.y), witness=witness,
                    consequence_equations=equations, equations=equations)

    def assert_exact_candidate(self, result, equations):
        self.assertEqual(result["status"], "ALGEBRAIC_CANDIDATE", result.get("reason"))
        candidate = result["candidate"]
        alpha = candidate["primitive_symbol"]
        modulus = sp.Poly(candidate["primitive_polynomial"], alpha, domain=sp.QQ)
        self.assertFalse(modulus.as_expr().has(sp.Float))
        for expression in candidate["squared_radii"]:
            self.assertFalse(expression.has(sp.Float))
        for equation in equations:
            expression = sp.cancel(equation.subs(candidate["radius_expressions"]))
            numerator, denominator = expression.as_numer_denom()
            # These rational-coefficient controls can be replayed directly
            # in the proposed primitive field, without floating evaluation.
            self.assertTrue(sp.Poly(numerator, alpha, domain=sp.QQ).rem(modulus).is_zero)
            self.assertFalse(sp.Poly(denominator, alpha, domain=sp.QQ).rem(modulus).is_zero)
        return candidate

    def test_homogeneous_rational_substitution_is_exact(self):
        x, y, z = sp.symbols("x y z")
        equation = x**4 + y*x**2 - 3*x + 2
        replacement = (y + z*z)/(y - z)
        cleared = _substitute(equation, {x: replacement})
        self.assertEqual(sp.cancel(cleared/(y-z)**4 - equation.subs(x, replacement)), 0)

    def test_irrational_radius_difference(self):
        x, y = self.x, self.y
        equations = [x+y-3, (y-x)**2-2]
        model = self.model(equations, {x: (3-2**0.5)/2, y: (3+2**0.5)/2})
        result = solve_model(model, max_seconds=2)
        candidate = self.assert_exact_candidate(result, equations)
        self.assertEqual(sp.degree(candidate["primitive_polynomial"]), 2)

    def test_exact_radical_base_coefficient(self):
        x, y = self.x, self.y
        result = solve_model(self.model([x-1, y-sp.sqrt(2)], {x: 1.0, y: 2**0.5}),
                             max_seconds=2)
        candidate = self.assert_exact_candidate(result, [x-1, y*y-2])
        self.assertGreater(float(candidate["algebraic_values"][y]), 0)

    def test_witness_does_not_supply_algebraic_coefficients(self):
        x, y = self.x, self.y
        equations = [x-1, y*y-2]
        with patch.object(sp, "nsimplify", side_effect=AssertionError("numerical recognition forbidden")):
            result = solve_model(self.model(equations, {x: 1.0+1e-9, y: 2**0.5+1e-9}),
                                 max_seconds=2)
            self.assert_exact_candidate(result, equations)
            bad = solve_model(self.model(equations, {x: 1.0, y: 1.4}), max_seconds=2)
        self.assertEqual(bad["status"], "UNRESOLVED")
        self.assertEqual(bad["candidates"], [])

    def test_tiny_budget_stays_unresolved(self):
        x, y = self.x, self.y
        result = solve_model(self.model([x+y-3, (y-x)**2-2],
                                        {x: (3-2**0.5)/2, y: (3+2**0.5)/2}),
                             max_seconds=0.000001)
        self.assertEqual(result["status"], "UNRESOLVED")
        self.assertTrue(result.get("timed_out"))
        self.assertEqual(result["candidates"], [])


if __name__ == "__main__":
    unittest.main()
