"""Regression coverage for exact roots over algebraic coefficient fields."""
import unittest
from unittest import mock

import sympy as sp

from reduce import _real_roots_near, _time_limit


class NumberFieldRootTests(unittest.TestCase):
    def setUp(self):
        self.x=sp.Symbol('x',real=True)
        # An accidental expensive compositum must not hang this small suite.
        self.budget=_time_limit(10)
        self.budget.__enter__()
        self.addCleanup(self.budget.__exit__,None,None,None)

    def exact_equal(self,left,right):
        self.assertEqual(sp.to_number_field(left-right).coeffs(),[0])

    def unique_root(self,polynomial,expected,*,approximation=None,tolerance=1e-5):
        guide=float(sp.N(expected,20)) if approximation is None else approximation
        roots=_real_roots_near(polynomial,self.x,guide,tolerance=tolerance)
        self.assertEqual(len(roots),1)
        self.exact_equal(roots[0],expected)
        expression=polynomial.as_expr() if isinstance(polynomial,sp.Poly) else polynomial
        self.exact_equal(expression.subs(self.x,roots[0]),0)
        return roots[0]

    def test_saved_j20_terminal_quadratic(self):
        x=self.x
        polynomial=x*x-sp.Rational(5,2)*x+sp.sqrt(5)/10+sp.Rational(17,16)
        root=self.unique_root(polynomial,sp.Rational(5,4)+sp.sqrt((5-sp.sqrt(5))/10))
        self.assertEqual(root.poly.degree(),4)

    def test_rational_path_never_calls_number_field_norm(self):
        x=self.x
        with mock.patch.object(sp.Poly,'norm',side_effect=AssertionError('rational path called norm')):
            roots=_real_roots_near(x*x-2,x,1.4142135623730951)
            self.assertEqual(len(roots),1)
            self.exact_equal(roots[0],sp.sqrt(2))
            self.assertEqual(_real_roots_near((x-3)**2,x,3.),[sp.Integer(3)])

    def test_mixed_rational_and_radical_coefficients(self):
        x=self.x
        polynomial=x*x+(1-2*sp.sqrt(2))*x-sp.sqrt(2)
        self.unique_root(polynomial,1+sp.sqrt(2))
        self.unique_root(polynomial,sp.sqrt(2)-2)
        self.unique_root((3+sp.sqrt(5))*(x-2-sp.sqrt(5)),2+sp.sqrt(5))

    def test_poly_input_and_negative_coefficient_embedding(self):
        x=self.x
        polynomial=sp.Poly(x*x+2*x+sp.Rational(1,2)-sp.sqrt(5)/10,x,extension=True)
        self.unique_root(polynomial,-1+sp.sqrt((5+sp.sqrt(5))/10))

    def test_cubic_coefficient_field(self):
        x=self.x
        cube=sp.real_root(2,3)
        root=self.unique_root(x*x-1-cube,sp.sqrt(1+cube))
        self.assertEqual(root.poly.degree(),6)

    def test_multiquadratic_coefficient_field(self):
        x=self.x
        expected=1+sp.sqrt(2)+sp.sqrt(3)
        root=self.unique_root(x-expected,expected)
        self.assertEqual(root.poly.degree(),4)

    def test_rational_root_and_repeated_field_norm(self):
        x=self.x
        polynomial=(x-1)*(x-sp.sqrt(2))
        self.assertEqual(self.unique_root(polynomial,sp.Integer(1)),1)
        self.unique_root((x-sp.sqrt(2))**2,sp.sqrt(2))

    def test_exact_embedding_rejects_arbitrarily_close_wrong_conjugate(self):
        x=self.x
        epsilon=sp.Rational(1,10**12)
        expected=1+epsilon*sp.sqrt(5)
        wrong=1-epsilon*sp.sqrt(5)
        polynomial=x-expected
        # Both roots lie in the requested window, and the wrong conjugate
        # has residual <1e-6; numerical residual filtering would accept it.
        self.assertLess(abs(float(polynomial.subs(x,wrong))),1e-6)
        self.unique_root(polynomial,expected,approximation=float(wrong))

    def test_unsupported_exact_embedding_remains_unresolved(self):
        for error in (sp.polys.polyerrors.IsomorphismFailed,
                      sp.polys.polyerrors.NotAlgebraic):
            with self.subTest(error=error.__name__):
                with mock.patch.object(sp,'to_number_field',side_effect=error('unsupported embedding')):
                    with self.assertRaisesRegex(ValueError,'coefficient embedding exactly: '+error.__name__):
                        _real_roots_near(self.x-sp.sqrt(2),self.x,1.4142135623730951)


if __name__=='__main__':
    unittest.main()
