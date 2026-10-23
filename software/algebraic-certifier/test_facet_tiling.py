"""Exact planar patch checks after face convexity has been established.

These isolate malformed arrangements which must not be accepted merely because
all their polygons lie in the same supporting plane.  Individual polygons in
the fixtures are convex; their metric regularity is an upstream obligation.
"""
import unittest

import sympy as sp

from verify import ExactField, ProofFailure
from verify_gram import _certify_facet_patch


class FacetTilingTests(unittest.TestCase):
    def certify(self, coordinates, faces, coplanar_vertices=None):
        points = [sp.Matrix([x, y, 0]) for x, y in coordinates]
        if coplanar_vertices is None:
            coplanar_vertices = list(range(len(points)))
        return _certify_facet_patch(list(range(len(faces))), faces, points,
                                   coplanar_vertices, ExactField())

    def test_two_triangles_tile_a_rhombus(self):
        # In a basis with mutual inner product 1/2, these are unit equilateral
        # triangles whose union is a rhombus with a shared short diagonal.
        coordinates = [(0, 0), (1, 0), (1, 1), (0, 1)]
        report = self.certify(coordinates, [[0, 1, 3], [1, 2, 3]])
        self.assertEqual(report['boundary'], [0, 1, 2, 3])
        self.assertEqual(report['internal_edges'], [[1, 3]])

    def test_independent_face_reversal_and_order_do_not_change_tiling(self):
        coordinates = [(0, 0), (1, 0), (1, 1), (0, 1)]
        for faces in ([[3, 1, 0], [1, 2, 3]],
                      [[3, 2, 1], [0, 1, 3]],
                      [[3, 2, 1], [3, 1, 0]]):
            with self.subTest(faces=faces):
                report = self.certify(coordinates, faces)
                self.assertEqual(set(report['boundary']), set(range(4)))
                self.assertEqual(report['internal_edges'], [[1, 3]])

    def test_three_convex_polygons_tile_a_strictly_convex_boundary(self):
        coordinates = [(0, 0), (2, 0), (3, 1), (1, 3), (0, 1)]
        report = self.certify(coordinates,
                              [[0, 1, 2], [0, 2, 3], [0, 3, 4]])
        self.assertEqual(report['boundary'], list(range(5)))
        self.assertEqual(report['internal_edges'], [[0, 2], [0, 3]])

    def test_same_side_overlapping_triangles_are_rejected(self):
        coordinates = [(0, 0), (2, 0), (2, 2), (0, 2)]
        with self.assertRaisesRegex(ProofFailure, 'overlap|inconsistent'):
            self.certify(coordinates, [[0, 1, 2], [0, 1, 3]])

    def test_duplicate_polygon_with_reversed_input_order_is_rejected(self):
        with self.assertRaisesRegex(ProofFailure, 'overlap|inconsistent'):
            self.certify([(0, 0), (1, 0), (0, 1)],
                         [[0, 1, 2], [2, 1, 0]])

    def test_missing_part_of_full_supporting_facet_is_rejected(self):
        coordinates = [(0, 0), (1, 0), (1, 1), (0, 1)]
        with self.assertRaisesRegex(ProofFailure, 'full supporting facet'):
            self.certify(coordinates, [[0, 1, 3]])

    def test_disconnected_boundary_cycles_are_rejected(self):
        coordinates = [(0, 0), (1, 0), (0, 1),
                       (3, 0), (4, 0), (3, 1)]
        with self.assertRaisesRegex(ProofFailure, 'multiple boundary'):
            self.certify(coordinates, [[0, 1, 2], [3, 4, 5]])

    def test_hole_is_rejected_despite_consistent_internal_edges(self):
        coordinates = [(-2, -2), (2, -2), (2, 2), (-2, 2),
                       (-1, -1), (1, -1), (1, 1), (-1, 1)]
        faces = [[0, 1, 5, 4], [1, 2, 6, 5],
                 [2, 3, 7, 6], [3, 0, 4, 7]]
        with self.assertRaisesRegex(ProofFailure, 'multiple boundary'):
            self.certify(coordinates, faces)

    def test_nonconvex_outer_boundary_is_rejected(self):
        coordinates = [(0, 0), (2, 0), (sp.Rational(1, 2),
                       sp.Rational(1, 2)), (0, 2)]
        with self.assertRaisesRegex(ProofFailure, 'not strictly convex'):
            self.certify(coordinates, [[0, 1, 2], [0, 2, 3]])

    def test_boundary_branches_at_a_shared_vertex(self):
        coordinates = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)]
        with self.assertRaisesRegex(ProofFailure, 'branches|repeats'):
            self.certify(coordinates, [[0, 1, 2], [0, 3, 4]])

    def test_collinear_boundary_vertex_is_rejected_by_curvature_contract(self):
        coordinates = [(0, 0), (1, 0), (2, 0), (0, 1)]
        with self.assertRaisesRegex(ProofFailure, 'not strictly convex'):
            self.certify(coordinates, [[0, 1, 3], [1, 2, 3]])


if __name__ == '__main__':
    unittest.main()
