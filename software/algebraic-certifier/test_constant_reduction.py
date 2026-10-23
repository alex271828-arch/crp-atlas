"""Exact constant closure and constant factors in radical reduction."""
import unittest

import sympy as sp

from reduce import (_groebner_solutions, _numerator, _radical_square_classes,
                    _time_limit, _triangular_solutions, _unique, solve_low_degree)


class ConstantReductionTests(unittest.TestCase):
    def setUp(self):
        self.budget = _time_limit(10)
        self.budget.__enter__()
        self.addCleanup(self.budget.__exit__, None, None, None)
        # These positive radicals are equal, but expansion does not prove it.
        self.zero = sp.sqrt(6 + 2*sp.sqrt(5)) - 1 - sp.sqrt(5)

    def assert_exact_equal(self, left, right):
        self.assertEqual(sp.to_number_field(left - right).coeffs(), [0])

    def test_algebraic_zero_constant_is_removed(self):
        self.assertNotEqual(self.zero.expand(), 0)
        self.assertEqual(_numerator(self.zero), 0)
        self.assertEqual(_unique([self.zero]), [])

    def test_algebraic_denominator_cancellation(self):
        root = sp.sqrt(5)
        residual = (15 + 7*root)/(4 + 4*root) - (5 + 2*root)/4
        self.assertNotEqual(residual.expand(), 0)
        self.assertEqual(_numerator(residual), 0)
        self.assertEqual(_groebner_solutions([residual], [], {}, []), [{}])

    def test_tiny_exact_nonzero_constant_is_retained(self):
        for exponent in (100, 300):
            tiny = sp.Rational(1, 10**exponent)
            for residual in (tiny, self.zero + tiny, self.zero - tiny):
                with self.subTest(residual=residual):
                    self.assertNotEqual(_numerator(residual), 0)
                    self.assertEqual(len(_unique([residual])), 1)

    def test_rootof_identity_closes_without_algebraic_assumption(self):
        x = sp.Symbol("x")
        root = sp.CRootOf(x**3 - x - 1, 0)
        identity = root**3 - root - 1
        self.assertIsNone(identity.is_algebraic)
        self.assertEqual(_numerator(identity), 0)
        self.assertEqual(_groebner_solutions([identity], [], {}, []), [{}])
        residual = identity + sp.Rational(1, 10**100)
        self.assertNotEqual(_numerator(residual), 0)
        self.assertEqual(_groebner_solutions([residual], [], {}, []), [])

    def test_no_variable_closure_accepts_algebraic_zero(self):
        for equations in ([], [sp.S.Zero], [self.zero], [sp.S.Zero, self.zero]):
            with self.subTest(equations=equations):
                self.assertEqual(_groebner_solutions(equations, [], {}, []), [{}])

    def test_no_variable_closure_rejects_tiny_exact_contradictions(self):
        tiny = sp.Rational(1, 10**100)
        for residual in (tiny, self.zero + tiny, self.zero - tiny):
            with self.subTest(residual=residual):
                self.assertEqual(_groebner_solutions([residual], [], {}, []), [])

    def test_triangular_no_variable_closure_checks_exact_constants(self):
        for equations in ([], [sp.S.Zero], [self.zero], [sp.S.Zero, self.zero]):
            with self.subTest(equations=equations):
                self.assertEqual(_triangular_solutions(equations, [], {}, []), [{}])
        tiny = sp.Rational(1, 10**100)
        for residual in (tiny, self.zero + tiny, self.zero - tiny):
            with self.subTest(residual=residual):
                self.assertEqual(_triangular_solutions([residual], [], {}, []), [])

    def test_linear_elimination_closes_an_algebraic_identity(self):
        x = sp.Symbol("x", real=True)
        expected = 1 + sp.sqrt(5)
        equations = [x - sp.sqrt(6 + 2*sp.sqrt(5)), x - expected]
        result = solve_low_degree(equations, [x], {x: float(expected)}, max_seconds=3)
        self.assertEqual(result["status"], "ALGEBRAIC_CANDIDATE", result.get("reason"))
        self.assertEqual(len(result["solutions"]), 1)
        self.assert_exact_equal(result["solutions"][0][x], expected)
        self.assertEqual(result["remaining_equations"], [])

    def test_linear_elimination_preserves_tiny_exact_contradiction(self):
        x = sp.Symbol("x", real=True)
        expected = 1 + sp.sqrt(5)
        equations = [x - sp.sqrt(6 + 2*sp.sqrt(5)),
                     x - expected + sp.Rational(1, 10**100)]
        result = solve_low_degree(equations, [x], {x: float(expected)}, max_seconds=3)
        self.assertEqual(result["status"], "UNRESOLVED")
        self.assertEqual(result["solutions"], [])
        self.assertTrue(result["remaining_equations"])

    def test_positive_algebraic_constant_radicals(self):
        y = sp.Symbol("y", positive=True)
        for radicand in ((21 + 9*sp.sqrt(5))/8, 2*sp.sqrt(5) + 6):
            with self.subTest(radicand=radicand):
                model = {"corner_types": [{"y": y, "S": radicand}]}
                replacements, equations, auxiliaries = _radical_square_classes(
                    model, {}, {y: float(sp.sqrt(radicand))}, [])
                self.assertEqual(set(replacements), {y})
                self.assert_exact_equal(replacements[y]**2, radicand)
                self.assertGreater(float(replacements[y]), 0)
                self.assertFalse(replacements[y].has(sp.Float))
                self.assertEqual(equations, [])
                self.assertEqual(auxiliaries, ())

    def test_negative_constant_radicands_do_not_supply_real_replacements(self):
        y = sp.Symbol("y", positive=True)
        for radicand in (-6 + 2*sp.sqrt(5), -12 + 4*sp.sqrt(5)):
            with self.subTest(radicand=radicand):
                model = {"corner_types": [{"y": y, "S": radicand}]}
                replacements, equations, auxiliaries = _radical_square_classes(
                    model, {}, {}, [])
                self.assertEqual(replacements, {})
                self.assertEqual(equations, [])
                self.assertEqual(auxiliaries, ())

    def test_constant_factors_keep_exact_positive_symbolic_radicals(self):
        x = sp.Symbol("x", real=True)
        y = sp.Symbol("y", positive=True)
        cases = ((x/8, 2),
                 ((x + 1)**2*(x + 2)/8, 2),
                 (x/(2*sp.sqrt(5) + 6), 2),
                 ((-6 + 2*sp.sqrt(5))*x, -1),
                 (x/(-12 + 4*sp.sqrt(5)), -1))
        for radicand, x_value in cases:
            with self.subTest(radicand=radicand):
                model = {"corner_types": [{"y": y, "S": radicand}]}
                witness = {x: float(x_value)}
                witness[y] = float(sp.sqrt(radicand.subs(x, x_value)))
                replacements, equations, auxiliaries = _radical_square_classes(
                    model, {}, witness, [])
                self.assertEqual(set(replacements), {y})
                self.assertEqual(len(auxiliaries), 1)
                self.assertEqual(len(equations), 1)
                residual = sp.cancel(replacements[y]**2 - radicand).as_numer_denom()[0]
                remainder = sp.rem(residual, equations[0], auxiliaries[0])
                self.assertEqual(sp.simplify(remainder), 0)
                self.assertGreater(float(replacements[y].subs(witness)), 0)
                self.assertFalse(replacements[y].has(sp.Float))


if __name__ == "__main__":
    unittest.main()
