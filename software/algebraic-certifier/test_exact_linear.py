"""Exact-field arithmetic must improve representation without changing equations."""
import unittest
from pathlib import Path
import subprocess
import sys

import sympy as sp
from sympy.matrices.exceptions import NonInvertibleMatrixError, NonSquareMatrixError
from sympy.polys.polyerrors import PolynomialError

from exact_linear import exact_inverse, normalize_polynomials


class PackageImportTests(unittest.TestCase):
    def test_gram_model_imports_helper_through_its_package(self):
        directory = Path(__file__).resolve().parent
        code = "\n".join([
            "import importlib, sys, types",
            "package = types.ModuleType('solver_package')",
            f"package.__path__ = [{str(directory)!r}]",
            "sys.modules[package.__name__] = package",
            "module = importlib.import_module('solver_package.gram_model')",
            "assert module.exact_inverse.__module__ == 'solver_package.exact_linear'",
        ])
        # Isolated mode excludes cwd/PYTHONPATH, so a broken absolute import
        # cannot accidentally find the script directory used by this test run.
        result = subprocess.run([sys.executable, "-I", "-c", code],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)


class ExactInverseTests(unittest.TestCase):
    def assert_inverse(self, matrix):
        inverse = exact_inverse(matrix)
        residual = matrix * inverse - sp.eye(matrix.rows)
        self.assertEqual(normalize_polynomials(list(residual), ()),
                         [sp.S.Zero] * len(residual))
        self.assertFalse(inverse.has(sp.Float))
        return inverse

    def test_rational_matrix(self):
        self.assertEqual(self.assert_inverse(sp.Matrix([[2, 1], [3, 2]])),
                         sp.Matrix([[2, -1], [-3, 2]]))

    def test_quadratic_coordinate_change_has_reduced_coefficients(self):
        q = sp.sqrt(5)
        matrix = sp.Matrix([[1, 1 + q, 3 - q],
                            [q - 1, 2, 1/(1 + q)],
                            [3, 2*q, 5 - q]])
        inverse = self.assert_inverse(matrix)
        for value in inverse:
            # The inverse uses the basis 1, sqrt(5), with rational coefficients;
            # no nested quotient containing a radical remains.
            self.assertTrue(sp.denom(value).is_Rational)
            self.assertLessEqual(sp.Poly(value.xreplace({q: sp.Symbol("q")}),
                                        sp.Symbol("q"), domain=sp.QQ).degree(), 1)

    def test_multiquadratic_matrix(self):
        self.assert_inverse(sp.Matrix([[sp.sqrt(2), 1 + sp.sqrt(3)],
                                       [1 - sp.sqrt(3), 2 + sp.sqrt(2)]]))

    def test_algebraically_singular_matrix(self):
        a = sp.sqrt(6 + 2*sp.sqrt(5))
        with self.assertRaises(NonInvertibleMatrixError):
            exact_inverse(sp.Matrix([[a, 2*a], [1 + sp.sqrt(5), 2 + 2*sp.sqrt(5)]]))

    def test_non_square_and_empty_matrices(self):
        with self.assertRaises(NonSquareMatrixError):
            exact_inverse(sp.Matrix([[1, 2]]))
        self.assertEqual(exact_inverse(sp.zeros(0)), sp.zeros(0))

    def test_rejects_inexact_or_nonconstant_coefficients(self):
        for value in (sp.Float("1.4142135623730951"), sp.Symbol("t"), sp.pi):
            with self.subTest(value=value), self.assertRaises(ValueError):
                exact_inverse(sp.Matrix([[value]]))


class NormalizePolynomialTests(unittest.TestCase):
    def setUp(self):
        self.x, self.y = sp.symbols("x y", real=True)
        self.q = sp.sqrt(5)
        self.zero = sp.sqrt(6 + 2*self.q) - 1 - self.q

    def test_simultaneous_affine_normalization_removes_algebraic_zero(self):
        x, y, q = self.x, self.y, self.q
        expressions = [(x + 2*y)/(1 + q),
                       (self.zero + 2)*x + self.zero*y,
                       self.zero]
        self.assertEqual(normalize_polynomials(expressions, (x, y)),
                         [x*(-sp.Rational(1, 4) + q/4)
                          + y*(-sp.Rational(1, 2) + q/2), 2*x, 0])

    def test_polynomial_powers_and_multiquadratic_coefficients(self):
        x, y = self.x, self.y
        values = normalize_polynomials(
            [(sp.sqrt(2) + sp.sqrt(3))*(x + y)**2,
             x/(sp.sqrt(2) + sp.sqrt(3))], (x, y))
        self.assertEqual(sp.Poly(values[0], x, y).total_degree(), 2)
        self.assertEqual(values[1], x*(-sp.sqrt(2) + sp.sqrt(3)))
        self.assertEqual(normalize_polynomials(
            [values[0] - (sp.sqrt(2) + sp.sqrt(3))*(x + y)**2], (x, y)), [0])

    def test_constants_and_empty_input(self):
        self.assertEqual(normalize_polynomials([self.zero, 1/(1 + self.q)], ()),
                         [0, (-1 + self.q)/4])
        self.assertEqual(normalize_polynomials([], (self.x,)), [])

    def test_tiny_exact_nonzero_coefficient_is_preserved(self):
        tiny = sp.Rational(1, 10**300)
        self.assertEqual(normalize_polynomials([(self.zero + tiny)*self.x],
                                              (self.x,)), [tiny*self.x])

    def test_rejects_unlisted_symbolic_variables(self):
        with self.assertRaises(ValueError):
            normalize_polynomials([self.x + self.y], (self.x,))

    def test_rejects_floats_transcendentals_and_nonpolynomials(self):
        for expression in (sp.Float("1.0")*self.x, sp.pi*self.x):
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                normalize_polynomials([expression], (self.x,))
        with self.assertRaises(PolynomialError):
            normalize_polynomials([1/self.x], (self.x,))

    def test_rejects_invalid_variable_lists(self):
        with self.assertRaises(TypeError):
            normalize_polynomials([self.x], (self.x + 1,))
        with self.assertRaises(ValueError):
            normalize_polynomials([self.x], (self.x, self.x))


if __name__ == "__main__":
    unittest.main()
