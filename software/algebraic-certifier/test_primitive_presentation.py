"""Preserve natural closure parameters when their field contains the constants."""
import unittest
from unittest import mock

import sympy as sp

import reduce
from reduce import ReductionTimeout, _primitive_candidate, _time_limit


class PrimitivePresentationTests(unittest.TestCase):
    def setUp(self):
        self.x = sp.Symbol("x")
        self.alpha = sp.Symbol("alpha", real=True)
        self.budget = _time_limit(20)
        self.budget.__enter__()
        self.addCleanup(self.budget.__exit__, None, None, None)

    def candidate(self, values):
        return _primitive_candidate(values, sp.symbols("s0:" + str(len(values))), {}, [])

    def normal(self, expression, polynomial):
        numerator, denominator = sp.cancel(expression).as_numer_denom()
        numerator = sp.Poly(numerator, self.alpha, domain=sp.QQ)
        denominator = sp.Poly(denominator, self.alpha, domain=sp.QQ)
        return (numerator * sp.invert(denominator, polynomial)).rem(polynomial).as_expr()

    def test_degree_eight_closure_keeps_its_radius_difference(self):
        x, alpha = self.x, self.alpha
        polynomial = 25*x**8 + 150*x**6 - 110*x**4 + 10*x**2 + 1
        root = sp.CRootOf(polynomial, 3)
        q = sp.sqrt(5)
        base = (5*root**2 + 20 + 8*q)/(10 + 2*q)
        result = self.candidate([base, base + root, base + (1 + q)*root/2])
        self.assertEqual(result["primitive_value"], root)
        self.assertEqual(result["primitive_polynomial"], polynomial.subs(x, alpha))
        p = sp.Poly(result["primitive_polynomial"], alpha)
        q_expression = (9 - 15*alpha**2 - 5*alpha**4)/(7*alpha**2 - 4)
        base_expression = (5*alpha**6 + 30*alpha**4 - 5*alpha**2 - 12)/(8*(7*alpha**2 - 4))
        expected = [base_expression, base_expression + alpha,
                    base_expression + (1 + q_expression)*alpha/2]
        self.assertEqual(result["squared_radii"], [self.normal(v, p) for v in expected])
        self.assertEqual(result["derivation"][-1]["operation"],
                         "primitive_root_coefficient_embedding")

    def test_coefficient_embedding_selects_the_correct_real_conjugate(self):
        x, alpha = self.x, self.alpha
        polynomial = x**4 - 10*x**2 + 1
        # These two positive roots give opposite polynomial expressions for
        # the same positive sqrt(2).  A field-isomorphism-only check is not
        # enough; the real coefficient embedding must also be respected.
        for index, sign in ((3, 1), (2, -1)):
            with self.subTest(index=index):
                root = sp.CRootOf(polynomial, index)
                result = self.candidate([root + sp.sqrt(2),
                                         (sp.sqrt(2) + 1)/(root + 2)])
                self.assertEqual(result["primitive_value"], root)
                p = sp.Poly(result["primitive_polynomial"], alpha)
                q = sign*(alpha**3 - 9*alpha)/2
                self.assertEqual(result["squared_radii"],
                                 [self.normal(alpha + q, p),
                                  self.normal((q + 1)/(alpha + 2), p)])

    def test_rational_coefficient_fast_path_is_unchanged(self):
        root = sp.CRootOf(self.x**3 - self.x - 1, 0)
        with mock.patch.object(reduce, "_root_field_expressions",
                               side_effect=AssertionError("optional embedding called")), \
                mock.patch.object(sp, "to_number_field",
                                  side_effect=AssertionError("compositum called")):
            result = self.candidate([root, 1/(root + 1), root**2 + 2])
        self.assertEqual(result["primitive_value"], root)
        p = sp.Poly(result["primitive_polynomial"], self.alpha)
        self.assertEqual(result["squared_radii"],
                         [self.normal(v, p) for v in
                          (self.alpha, 1/(self.alpha + 1), self.alpha**2 + 2)])
        self.assertEqual(result["derivation"], [])

    def test_radius_in_proper_subfield_keeps_its_smaller_degree(self):
        root = sp.CRootOf(self.x**4 - 10*self.x**2 + 1, 3)
        result = self.candidate([root**2 + sp.sqrt(6)])
        self.assertNotEqual(result["primitive_value"], root)
        self.assertEqual(result["primitive_polynomial"],
                         self.alpha**2 - 10*self.alpha - 29)
        self.assertEqual(result["squared_radii"], [self.alpha])
        self.assertEqual(result["derivation"], [])

    def test_radii_collectively_generate_full_field(self):
        root = sp.CRootOf(self.x**4 - 10*self.x**2 + 1, 3)
        # Each radius lies in a quadratic subfield; together they generate
        # the full quartic field.  Testing each radius separately would miss
        # this valid preservation of the closure parameter.
        values = [(root**3 - 9*root)/2 + sp.sqrt(2),
                  (11*root - root**3)/2 + sp.sqrt(3)]
        result = self.candidate(values)
        self.assertEqual(result["primitive_value"], root)
        self.assertEqual(result["squared_radii"],
                         [self.alpha**3 - 9*self.alpha,
                          11*self.alpha - self.alpha**3])
        self.assertEqual(result["derivation"][-1]["operation"],
                         "primitive_root_coefficient_embedding")

    def test_coefficient_outside_root_field_uses_compositum(self):
        root = sp.CRootOf(self.x**2 - 2, 1)
        result = self.candidate([root + sp.sqrt(3)])
        self.assertNotEqual(result["primitive_value"], root)
        self.assertEqual(sp.Poly(result["primitive_polynomial"], self.alpha).degree(), 4)
        self.assertEqual(result["squared_radii"], [self.alpha])
        self.assertEqual(result["derivation"], [])

    def test_optional_embedding_timeout_falls_back_to_compositum(self):
        root = sp.CRootOf(self.x**2 - 2, 1)
        with mock.patch.object(reduce, "_root_field_expressions",
                               side_effect=ReductionTimeout("optional attempt expired")) as attempt, \
                mock.patch.object(reduce, "_time_limit", wraps=_time_limit) as timer:
            result = self.candidate([root + sp.sqrt(3)])
        attempt.assert_called_once()
        timer.assert_called_once_with(2.0)
        self.assertEqual(sp.Poly(result["primitive_polynomial"], self.alpha).degree(), 4)
        self.assertEqual(result["squared_radii"], [self.alpha])

    def test_unsupported_optional_embedding_falls_back_to_compositum(self):
        root = sp.CRootOf(self.x**2 - 2, 1)
        with mock.patch.object(reduce, "_root_field_expressions",
                               side_effect=NotImplementedError("unsupported embedding")):
            result = self.candidate([root + sp.sqrt(3)])
        self.assertEqual(sp.Poly(result["primitive_polynomial"], self.alpha).degree(), 4)

    def test_full_field_proof_timeout_falls_back_to_compositum(self):
        root = sp.CRootOf(self.x**2 - 2, 1)
        with mock.patch.object(reduce, "_radii_generate_root_field",
                               side_effect=ReductionTimeout("field proof expired")) as proof, \
                mock.patch.object(reduce, "_time_limit", wraps=_time_limit) as timer:
            result = self.candidate([root + sp.sqrt(2)])
        proof.assert_called_once()
        timer.assert_called_once_with(2.0)
        self.assertNotEqual(result["primitive_value"], root)
        self.assertEqual(result["primitive_polynomial"], self.alpha**2 - 8)
        self.assertEqual(result["squared_radii"], [self.alpha])


if __name__ == "__main__":
    unittest.main()
