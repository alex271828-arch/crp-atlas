"""Projection completion must retain algebraic regular-face constants."""
import unittest
import sympy as sp
from gram_reduction import eliminate_affine_equations

class AffineProjectionTests(unittest.TestCase):
    def test_explicit_base_field_coordinate_needs_no_root_lifting(self):
        x=sp.Symbol('x',real=True)
        value=(1+sp.sqrt(5))/4
        equations,variables,assignment=eliminate_affine_equations([x-value],[x],[])
        self.assertEqual(equations,[])
        self.assertEqual(variables,[])
        self.assertEqual(assignment,{x:value})

    def test_affine_completion_preserves_remaining_nonlinear_geometry(self):
        x,y=sp.symbols('x y',real=True)
        value=(1+sp.sqrt(5))/4
        equations,variables,assignment=eliminate_affine_equations([x-value,y*y-x],[x,y],[])
        self.assertEqual(variables,[y])
        self.assertEqual(assignment,{x:value})
        self.assertEqual(len(equations),1)
        self.assertEqual(sp.expand(equations[0]-(y*y-value)),0)

    def test_coupled_affine_system_over_algebraic_field(self):
        x,y=sp.symbols('x y',real=True)
        equations,variables,assignment=eliminate_affine_equations(
            [sp.sqrt(5)*x+y-2,x-y-1],[x,y],[])
        self.assertEqual(equations,[])
        self.assertEqual(variables,[])
        self.assertEqual(sp.simplify(sp.sqrt(5)*assignment[x]+assignment[y]-2),0)
        self.assertEqual(sp.simplify(assignment[x]-assignment[y]-1),0)

if __name__=='__main__':unittest.main()
