"""Exact metric-patch tests, including false-square and relabeling guards."""

import copy
import unittest

import sympy as sp

try:
    from .patches import (combinatorial_automorphisms, enrich_system,
                          infer_reflection_squares, metric_patch_equation,
                          squared_polygon_distance)
except ImportError:
    from patches import (combinatorial_automorphisms, enrich_system,
                         infer_reflection_squares, metric_patch_equation,
                         squared_polygon_distance)


TETRA = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]]
PYRAMID = [[0, 1, 2, 3], [4, 1, 0], [4, 2, 1], [4, 3, 2], [4, 0, 3]]
BIPYRAMID = [[0, 1, 2], [0, 2, 3], [0, 3, 1], [4, 2, 1], [4, 3, 2], [4, 1, 3]]
J87 = [[6, 5, 10], [0, 2, 3], [5, 1, 7], [8, 6, 10], [9, 3, 8],
       [0, 4, 1], [0, 3, 4], [9, 8, 10], [3, 2, 8], [7, 1, 4],
       [7, 4, 9], [10, 7, 9], [8, 2, 6], [5, 7, 10], [0, 1, 2],
       [2, 1, 5, 6], [4, 3, 9]]
J87_GROUPS = [[1, 2], [0, 3, 4, 5, 6], [7, 8], [9, 10]]


def model(faces, groups=None):
    vertices = sorted(set(sum(faces, [])))
    groups = [[v] for v in vertices] if groups is None else groups
    return {"faces": copy.deepcopy(faces), "groups": copy.deepcopy(groups),
            "radii_symbols": list(sp.symbols(f"s0:{len(groups)}"))}


class PatchTests(unittest.TestCase):
    def test_exact_polygon_distances_have_no_float_coefficients(self):
        self.assertEqual(squared_polygon_distance(3, 1), 1)
        self.assertEqual(squared_polygon_distance(4, 2), 2)
        self.assertEqual(squared_polygon_distance(5, 2), (3 + sp.sqrt(5)) / 2)
        self.assertEqual(squared_polygon_distance(6, 3), 4)
        for n in (3, 4, 5, 6, 8, 10):
            for k in range(n):
                value = squared_polygon_distance(n, k)
                self.assertFalse(value.atoms(sp.Float))
                self.assertEqual(value, squared_polygon_distance(n, n - k))

    def test_tetrahedral_quadric_and_complete_automorphism_group(self):
        m = enrich_system(model(TETRA, [list(range(4))]))
        r = m["radii_symbols"][0]
        self.assertIn(r - sp.Rational(3, 8), m["patch_equations"])
        self.assertEqual(len(m["patch_metadata"]["automorphism_enumeration"]["permutations"]), 24)
        self.assertTrue(m["patch_metadata"]["automorphism_enumeration"]["complete"])
        self.assertEqual(m["reflections"], [])  # Tetrahedron mirrors fix only two vertices.
        self.assertFalse(any(expression.atoms(sp.Float) for expression in m["patch_equations"]))

    def test_square_pyramid_gives_linear_extra_point_constraint(self):
        m = enrich_system(model(PYRAMID))
        s0, s1, s2, s3, s4 = m["radii_symbols"]
        expected = s0 - s1 + s2 - s3
        self.assertIn(expected, m["patch_equations"])
        substitution = dict(zip(m["radii_symbols"], [sp.Rational(5, 8)] * 4 + [sp.Rational(1, 8)]))
        self.assertTrue(all(sp.simplify(e.subs(substitution)) == 0 for e in m["patch_equations"]))
        self.assertTrue(any(source["kind"] == "additional_point_affine_radius_relation"
                            for source in m["patch_sources"]))
        self.assertTrue(m["square_pyramid_patches"])

    def test_coplanar_four_points_do_not_supply_a_spatial_metric_patch(self):
        known = {(0, 1): 1, (1, 2): 1, (2, 3): 1, (0, 3): 1,
                 (0, 2): 2, (1, 3): 2}
        radii = dict(enumerate(sp.symbols("r0:4")))
        self.assertIsNone(metric_patch_equation((0, 1, 2, 3), known, radii))

    def test_integer_distance_input_remains_exact_and_relabeling_is_invariant(self):
        known = {(i, j): 1 for i in range(4) for j in range(i + 1, 4)}
        patch = metric_patch_equation((0, 1, 2, 3), known, dict(enumerate(sp.symbols("r0:4"))))
        self.assertIsNotNone(patch)
        self.assertFalse(patch["equation"].atoms(sp.Float))
        for faces in (TETRA, PYRAMID):
            original = model(faces)
            n = len(original["groups"])
            relabeled_faces = [[n - 1 - v for v in (face[1:] + face[:1])[::-1]]
                               for face in faces[::-1]]
            relabeled_groups = [[n - 1 - v for v in group] for group in original["groups"]]
            first = enrich_system(original)
            second = enrich_system(model(relabeled_faces, relabeled_groups))
            self.assertEqual(set(first["patch_equations"]), set(second["patch_equations"]))

    def test_bare_unit_four_cycle_does_not_become_a_square(self):
        m = enrich_system(model(BIPYRAMID))
        # 0-1-4-2-0 is a unit four-cycle, but its diagonal 1-2 has length one.
        self.assertEqual(m["known_distances"][(1, 2)], 1)
        self.assertNotIn((0, 4), m["known_distances"])
        self.assertEqual(m["reflection_squares"], [])
        self.assertTrue(m["reflections"])

    def test_unequal_transposed_edges_and_conflicting_diagonals_are_rejected(self):
        mirror = {"id": 0, "swapped_pairs": [[0, 1], [2, 3]]}
        known = {(0, 1): 1, (2, 3): 2, (0, 2): 1, (1, 3): 1}
        self.assertEqual(infer_reflection_squares(known, mirror), [])
        known[(2, 3)] = 1
        self.assertEqual(len(infer_reflection_squares(known, mirror)), 1)
        known[(0, 3)] = 1
        self.assertEqual(infer_reflection_squares(known, mirror), [])

    def test_j87_hidden_square_quadric_and_relabeling(self):
        permutation = [7, 2, 10, 0, 5, 8, 3, 9, 1, 6, 4]
        for relabeled in (False, True):
            faces, groups = J87, J87_GROUPS
            square = {1, 2, 3, 4}
            fixed = {0, 9, 10}
            if relabeled:
                faces = [[permutation[v] for v in face[::-1]] for face in J87[::-1]]
                groups = [[permutation[v] for v in group] for group in J87_GROUPS]
                square, fixed = {permutation[v] for v in square}, {permutation[v] for v in fixed}
            with self.subTest(relabeled=relabeled):
                original = copy.deepcopy(faces)
                m = enrich_system(model(faces, groups), positions="ignored, not geometry")
                a, b, _, _ = m["radii_symbols"]
                expected = sp.Poly(3 * (a - b)**2 - 4 * (a + b) + 4, a, b).monic().as_expr()
                self.assertIn(sp.expand(expected), m["patch_equations"])
                self.assertTrue(any(set(s["vertices"]) == square for s in m["reflection_squares"]))
                self.assertTrue(any(set(r["fixed_vertices"]) == fixed for r in m["reflections"]))
                self.assertEqual(m["faces"], original)
                self.assertFalse(m["patch_metadata"]["numerical_positions_used"])
                self.assertEqual(m["patch_metadata"]["distance_conflicts"], [])

    def test_open_or_zero_curvature_complex_does_not_supply_a_mirror_proof(self):
        subdivided_prism = ([list(range(6, 12))]
                            + [[i, (i + 1) % 6, (i + 1) % 6 + 6, i + 6] for i in range(6)]
                            + [[12, i, (i + 1) % 6] for i in range(6)])
        for faces in ([list(range(4))], [], subdivided_prism):
            result = combinatorial_automorphisms(faces)
            self.assertEqual(result["permutations"], [])
            self.assertFalse(result["complete"])

    def test_bounds_are_explicit_and_do_not_claim_complete_enumeration(self):
        result = combinatorial_automorphisms(TETRA, max_automorphisms=2)
        self.assertEqual(len(result["permutations"]), 2)
        self.assertFalse(result["complete"])
        m = enrich_system(model(PYRAMID), max_metric_patches=1)
        self.assertEqual(len(m["metric_patches"]), 1)
        self.assertTrue(m["patch_metadata"]["metric_patches_truncated"])


if __name__ == "__main__":
    unittest.main()
