"""Exact construction tests, including degeneracy and unselected branches."""
import unittest

import sympy as sp

from distance_kernel import (ConstructionLimit, complete_gram_entry,
                             distance_intersection, gram_extension,
                             gram_rank_relations, quadratic_branches,
                             reflect_across_anchor_hyperplane)


class DistanceKernelTests(unittest.TestCase):
    def assertZero(self, value):
        self.assertEqual(sp.cancel(value), 0)

    def test_three_unit_spheres_keep_both_tetrahedral_apices(self):
        anchors = [[0, 0, 0], [1, 0, 0], [sp.Rational(1, 2), sp.sqrt(3)/2, 0]]
        result = distance_intersection(anchors, [1, 1, 1])
        self.assertEqual(result["rank"], 2)
        self.assertEqual(len(result["free_parameters"]), 1)
        branches = quadratic_branches(result)["branches"]
        self.assertEqual(len(branches), 2)
        self.assertEqual({branch["point"][2] for branch in branches},
                         {-sp.sqrt(6)/3, sp.sqrt(6)/3})
        for branch in branches:
            for anchor in anchors:
                delta = branch["point"]-sp.Matrix(anchor)
                self.assertZero(delta.dot(delta)-1)

    def test_one_sphere_retains_a_surface(self):
        result = distance_intersection([[1, 2, 3]], [4])
        self.assertEqual(result["rank"], 0)
        self.assertEqual(len(result["free_parameters"]), 3)
        self.assertEqual(sp.Poly(result["residual_norm"], *result["free_parameters"]).total_degree(), 2)
        with self.assertRaises(ValueError):
            quadratic_branches(result)

    def test_overdetermined_inconsistency_is_retained(self):
        result = distance_intersection([[0, 0], [1, 0], [0, 1], [1, 1]], [1, 1, 1, 2])
        self.assertEqual(result["rank"], 2)
        self.assertEqual(len(result["free_parameters"]), 0)
        self.assertTrue(result["affine_residuals"])
        self.assertTrue(any(value.is_number and value != 0
                            for value in result["affine_residuals"]))

    def test_symbolic_pivot_keeps_its_nonzero_obligation(self):
        a = sp.Symbol("a", real=True)
        result = distance_intersection([[0, 0], [a, 0]], [1, 1], witness={a: 2.0})
        self.assertEqual(result["rank"], 1)
        self.assertTrue(any(condition.has(a) for condition in result["nonzero_conditions"]))
        self.assertTrue(result["provenance"]["rank_is_conditional_on_nonzero_pivot"])
        self.assertFalse(any(expression.has(sp.Float) for expression in result["affine_point"]))

    def test_rank_search_limit_does_not_misreport_a_lower_rank(self):
        with self.assertRaises(ConstructionLimit):
            distance_intersection([[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                                  [1, 1, 1], max_minors=1)

    def test_known_point_reflection_is_rational_and_preserves_distances(self):
        H = sp.Symbol("H", positive=True)
        metric = sp.diag(1, 1, H)
        anchors = [[0, 0, 0], [1, 0, 0], [0, 1, 1]]
        point = sp.Matrix([2, 3, 4])
        reflected = reflect_across_anchor_hyperplane(point, anchors, metric=metric)
        other = reflected["reflected_point"]
        self.assertEqual(reflected["known_point_equations"], [])
        self.assertFalse(any(e.has(sp.sqrt(H)) for e in other))
        self.assertNotEqual(other, point)
        for anchor in anchors:
            d0 = point-sp.Matrix(anchor)
            d1 = other-sp.Matrix(anchor)
            self.assertZero((d0.T*metric*d0)[0]-(d1.T*metric*d1)[0])
        twice = reflect_across_anchor_hyperplane(other, anchors, metric=metric)["reflected_point"]
        for value in twice-point:
            self.assertZero(value)

    def test_invalid_known_point_is_not_silently_accepted(self):
        result = distance_intersection([[0, 0], [1, 0]], [1, 1])
        branches = quadratic_branches(result, known_point=[sp.Rational(1, 2), 9])
        self.assertTrue(branches["known_point_equations"])

    def test_float_coefficients_are_rejected(self):
        with self.assertRaises(ValueError):
            distance_intersection([[0.1, 0, 0]], [1])

    def test_gram_schur_complement_recovers_orthogonal_height(self):
        x, y, z = sp.symbols("x y z", real=True)
        result = gram_extension(sp.eye(2), [x, y], x*x+y*y+z*z)
        self.assertZero(result["schur_residual"]-z*z)
        self.assertEqual(result["equations"], [])
        full = gram_extension(sp.eye(3), [x, y, z], x*x+y*y+z*z+1)
        self.assertEqual(full["equations"], [1])

    def test_gram_completion_keeps_the_two_orientation_choices(self):
        g = sp.Symbol("g", real=True)
        result = complete_gram_entry(sp.eye(2), [1, 2], [3, 4], 9, 50, unknown=g)
        self.assertZero(result["equation"]-((g-11)**2-100))
        self.assertEqual(set(result["branches"]), {1, 21})

    def test_rank_two_relations_hold_for_a_general_planar_gram_matrix(self):
        t = sp.Symbol("t", real=True)
        vectors = sp.Matrix([[1, 0, 2, t], [0, 1, 3, 4]])
        G = vectors.T*vectors
        result = gram_rank_relations(G, 2, witness={t: 1.5})
        self.assertEqual(result["constraints"], [])
        self.assertEqual(len(result["anchor_indices"]), 2)
        perturbed = G.copy()
        perturbed[3, 3] += 1
        result2 = gram_rank_relations(perturbed, 2, witness={t: 1.5})
        self.assertTrue(any(e != 0 for e in result2["constraints"]))

    def test_rank_one_conditions_are_quadratic_and_fraction_free(self):
        a, b, c, d, e, f = sp.symbols("a b c d e f", real=True)
        G = sp.Matrix([[a, b, c], [b, d, e], [c, e, f]])
        result = gram_rank_relations(G, 1, witness={a: 1, d: 0, f: 0})
        self.assertEqual(result["anchor_indices"], (0,))
        self.assertEqual(set(result["constraints"]), {a*d-b*b, a*e-b*c, a*f-c*c})
        self.assertEqual(result["positive_conditions"], [a])
        self.assertEqual(result["nonzero_conditions"], [a])

    def test_numeric_gram_guide_only_orders_exact_anchors(self):
        a, b, c = sp.symbols("a b c", real=True)
        G = sp.Matrix([[a, b], [b, c]])
        result = gram_rank_relations(G, 1, witness=[[0.0, 0.0], [0.0, 4.0]])
        self.assertEqual(result["anchor_indices"], (1,))
        self.assertEqual(result["constraints"], [a*c-b*b])
        self.assertFalse(any(e.has(sp.Float) for e in result["constraints"]))

    def test_rank_zero_has_no_nonzero_anchor_assumption(self):
        a, b, c = sp.symbols("a b c")
        result = gram_rank_relations([[a, b], [b, c]], 0)
        self.assertEqual(result["constraints"], [a, b, c])
        self.assertEqual(result["nonzero_conditions"], [])


if __name__ == "__main__":
    unittest.main()
