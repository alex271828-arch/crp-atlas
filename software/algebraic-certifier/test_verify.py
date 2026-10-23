"""Independent exact positive and negative tests of geometric certification."""
import unittest
from fractions import Fraction

import sympy as sp

try:
    from .verify import verify_solution, ExactField, RadicalAlgebra, Interval, ProofFailure, PI
except ImportError:
    from verify import verify_solution, ExactField, RadicalAlgebra, Interval, ProofFailure, PI


TETRA = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]]
CUBE = [[0, 1, 2, 3], [4, 7, 6, 5], [0, 4, 5, 1],
        [1, 5, 6, 2], [2, 6, 7, 3], [3, 7, 4, 0]]
BODY87 = [[12, 4, 13], [8, 7, 11], [11, 10, 6, 12], [9, 1, 10],
          [8, 11, 13], [0, 2, 3], [12, 6, 5], [3, 2, 7, 8],
          [9, 2, 1], [0, 3, 4, 5], [11, 12, 13], [7, 9, 10],
          [4, 3, 8], [0, 5, 6, 1], [7, 2, 9], [11, 7, 10],
          [13, 4, 8], [10, 1, 6], [0, 1, 2], [12, 5, 4]]
GROUPS87 = [[v for v in range(14) if v not in (9, 13)], [9, 13]]
BODY140 = [[18, 5, 19], [19, 4, 11, 20, 21], [22, 17, 16], [18, 19, 24],
           [24, 19, 21], [8, 7, 17], [15, 14, 18], [26, 20, 10, 29, 30],
           [22, 16, 25, 27, 28], [9, 2, 1, 12, 13], [23, 28, 31], [30, 29, 31],
           [13, 12, 23], [28, 27, 30], [4, 3, 11], [12, 8, 17, 22, 23],
           [15, 18, 24], [16, 15, 25], [0, 1, 2], [12, 1, 8], [14, 6, 5],
           [3, 2, 9, 10, 11], [19, 5, 4], [0, 6, 7, 8, 1], [18, 14, 5],
           [20, 11, 10], [0, 2, 3], [28, 30, 31], [0, 3, 4, 5, 6],
           [7, 14, 15, 16, 17], [29, 10, 9], [25, 24, 21, 26, 27],
           [13, 23, 31], [23, 22, 28], [7, 6, 14], [26, 21, 20],
           [29, 9, 13], [27, 26, 30], [25, 15, 24], [29, 13, 31]]


def prism_faces(n):
    return [list(range(n)), list(range(n, 2 * n))] + [
        [j, (j + 1) % n, (j + 1) % n + n, j + n] for j in range(n)]


class VerificationTests(unittest.TestCase):
    def assertCertified(self, report):
        self.assertTrue(report["certified"], report)
        self.assertEqual(report["status"], "CERTIFIED")

    def assertRejected(self, report):
        self.assertFalse(report["certified"], report)
        self.assertNotEqual(report["status"], "CERTIFIED")
        self.assertIn("reason", report)

    def test_exact_tetrahedron_passes_and_wrong_radii_fail(self):
        self.assertCertified(verify_solution(TETRA, [list(range(4))], [sp.Rational(3, 8)]))
        for radius in [sp.Rational(2, 5), sp.Rational(1, 4), 0, -1]:
            with self.subTest(radius=radius):
                self.assertRejected(verify_solution(TETRA, [list(range(4))], [radius]))

    def test_face_flatness_perturbation_fails(self):
        radii = [sp.Rational(3, 4)] * 8
        radii[0] += sp.Rational(1, 100)
        report = verify_solution(CUBE, [[v] for v in range(8)], radii)
        self.assertRejected(report)
        self.assertIn("flatness", report["reason"])

    def test_body87_coplanar_hinges_are_accepted(self):
        report = verify_solution(BODY87, GROUPS87, [1, 2])
        self.assertCertified(report)
        expected = {(4, 12), (8, 11), (1, 10), (1, 2), (7, 10), (4, 8)}
        self.assertEqual({tuple(edge) for edge in report["coplanar_edges"]}, expected)
        for cap_radius in [0, sp.Rational(1, 2), 1]:
            self.assertRejected(verify_solution(BODY87, GROUPS87, [1, cap_radius]))

    def test_positive_inward_cap_branch_rejected_by_convexity(self):
        a = (3 + sp.sqrt(5)) / 2
        groups = [[v for v in range(32) if v not in (18, 31)], [18, 31]]
        report = verify_solution(BODY140, groups, [a, (a + 1) / 5])
        self.assertRejected(report)
        self.assertIn("nonconvex hinge", report["reason"])

    def test_primitive_real_embedding_and_interval_are_respected(self):
        alpha = sp.Symbol("alpha")
        # The same symbolic expression is 3/8 in one real embedding and
        # 131/8 in the other. Both are positive; fan closure distinguishes them.
        radius = sp.Rational(3, 8) + (sp.sqrt(2) * alpha - 2) ** 2
        parameters = dict(primitive_symbol=alpha, primitive_polynomial=alpha ** 2 - 2)
        self.assertCertified(verify_solution(TETRA, [list(range(4))], [radius],
                                            isolating_interval=(sp.Rational(7, 5), sp.Rational(3, 2)), **parameters))
        self.assertRejected(verify_solution(TETRA, [list(range(4))], [radius],
                                           isolating_interval=(-sp.Rational(3, 2), -sp.Rational(7, 5)), **parameters))
        for interval in [(sp.Rational(1, 10), sp.Rational(1, 5)), (-2, 2), (2, 1)]:
            self.assertRejected(verify_solution(TETRA, [list(range(4))], [radius],
                                               isolating_interval=interval, **parameters))
        self.assertRejected(verify_solution(TETRA, [list(range(4))], [1 / (alpha ** 2 - 2)],
                                           isolating_interval=(1, 2), **parameters))

    def test_malformed_topology_and_groups_are_rejected(self):
        invalid = [(TETRA[:-1], [list(range(4))]),
                   (TETRA, [[0, 1, 2], [2, 3]]),
                   ([[0, 1, 1], *TETRA[1:]], [list(range(4))]),
                   (TETRA + [[v + 4 for v in f] for f in TETRA], [list(range(8))])]
        for faces, groups in invalid:
            self.assertRejected(verify_solution(faces, groups, [sp.Rational(3, 8)] * len(groups)))
        self.assertRejected(verify_solution([], [], []))

    def test_prisms_cover_all_regular_face_coefficient_fields(self):
        alpha = sp.Symbol("alpha")
        for n in (5, 6, 8, 10):
            with self.subTest(n=n):
                radius = sp.simplify(1 / (4 * sp.sin(sp.pi / n) ** 2) + sp.Rational(1, 4))
                polynomial = sp.minpoly(radius, alpha)
                low = sp.floor(10 * radius) / 10
                # Rational brackets select the intended real embedding;
                # root isolation is checked again inside the verifier.
                report = verify_solution(prism_faces(n), [list(range(2 * n))], [alpha],
                                         primitive_symbol=alpha, primitive_polynomial=polynomial,
                                         isolating_interval=(low, low + sp.Rational(1, 10)))
                self.assertCertified(report)
                self.assertEqual(report["checks"]["winding_number_one"], 2 * n)


class ExactArithmeticTests(unittest.TestCase):
    def test_rational_intervals_are_outward(self):
        value = Interval(Fraction(2), Fraction(3)).sqrt(digits=20)
        self.assertLessEqual(value.lo ** 2, 2)
        self.assertGreaterEqual(value.hi ** 2, 3)
        self.assertEqual((Interval(-2, 1) ** 2).lo, 0)
        self.assertEqual((Interval(-2, 1) ** 2).hi, 4)
        with self.assertRaises(ProofFailure):
            Interval(-1, 1).reciprocal()
        self.assertGreater(PI.lo, Fraction(314159265358979323846, 10 ** 20))
        self.assertLess(PI.hi, Fraction(314159265358979323847, 10 ** 20))

    def test_radical_square_classes_and_positive_roots(self):
        field = ExactField()
        algebra = RadicalAlgebra(field)
        root2, root8 = algebra.sqrt(2), algebra.sqrt(8)
        self.assertEqual(algebra.multiply(root2, root8), {0: sp.Integer(4)})
        root3, root6 = algebra.sqrt(3), algebra.sqrt(6)
        self.assertEqual(algebra.multiply(root2, root3), root6)
        self.assertEqual(len(algebra.generators), 2)
        with self.assertRaises(ProofFailure):
            algebra.sqrt(-1)
        alpha = sp.Symbol("alpha")
        field = ExactField(alpha, alpha ** 2 - 2, (1, 2))
        root = field.square_root_in_field((alpha - 2) ** 2)
        self.assertTrue(field.zero(root - (2 - alpha)))
        self.assertEqual(field.sign(root), 1)


if __name__ == "__main__":
    unittest.main()
