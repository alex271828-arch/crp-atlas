"""Small exact checks of the general encoder, independent of solver coverage."""
import unittest

import sympy as sp

try:
    from .encoder import build_system, _relation_consequences, squared_short_diagonal
except ImportError:
    from encoder import build_system, _relation_consequences, squared_short_diagonal


TETRA = [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)]
CUBE = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1),
        (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]
J88_FACES = [[7, 1, 8], [9, 3, 7], [4, 9, 10], [5, 4, 11],
             [11, 6, 5], [0, 2, 3], [0, 3, 4, 5], [3, 2, 7],
             [10, 8, 6], [11, 10, 6], [8, 1, 6], [9, 8, 10],
             [0, 5, 6, 1], [11, 4, 10], [0, 1, 2], [4, 3, 9],
             [9, 7, 8], [7, 2, 1]]
J88_GROUPS = [[1, 3, 4, 6], [8, 9], [0, 5, 7, 10], [2, 11]]
J88_RADII = [.6171077448892142, .8043821816376674,
             1.0968931800780984, 1.7353793334283981]


class EncoderTests(unittest.TestCase):
    def test_tetrahedron_exact_closure_and_radius(self):
        model = build_system(TETRA, [[0, 1, 2, 3]], [3 / 8])
        s, = model["radius_variables"]
        self.assertEqual(model["consequence_equations"], [8 * s - 3])
        assignment = {s: sp.Rational(3, 8)}
        for corner in model["corner_types"]:
            assignment[corner["y"]] = sp.sqrt(corner["S"].subs(assignment))
        self.assertTrue(all(sp.simplify(e.subs(assignment)) == 0 for e in model["equations"]))
        self.assertNotEqual(model["consequence_equations"][0].subs(s, sp.Rational(2, 5)), 0)

    def test_square_flatness_row_space_and_exact_cube_closure(self):
        model = build_system(CUBE, [[v] for v in range(8)], [.75] * 8)
        s = model["radius_variables"]
        rows, _ = sp.linear_eq_to_matrix(model["flatness_equations"], s)
        self.assertEqual(rows.rank(), 4)
        expected, _ = sp.linear_eq_to_matrix([s[0] - s[1] + s[2] - s[3]], s)
        self.assertEqual(rows.col_join(expected).rank(), rows.rank())
        assignment = {x: sp.Rational(3, 4) for x in s}
        for corner in model["corner_types"]:
            assignment[corner["y"]] = sp.sqrt(corner["S"].subs(assignment))
        self.assertTrue(all(sp.simplify(e.subs(assignment)) == 0 for e in model["equations"]))
        perturbed = dict(assignment, **{})
        perturbed[s[0]] += sp.Rational(1, 100)
        self.assertTrue(any(e.subs(perturbed) != 0 for e in model["flatness_equations"]))

    def test_regular_face_constants_are_exact(self):
        for n, value in [(3, 1), (4, 2), (5, (3 + sp.sqrt(5)) / 2),
                         (6, 3), (8, 2 + sp.sqrt(2)), (10, (5 + sp.sqrt(5)) / 2)]:
            self.assertEqual(sp.simplify(squared_short_diagonal(n) - value), 0)

    def test_face_order_orientation_and_vertex_names_do_not_matter(self):
        groups = [[0, 2, 4, 6], [1, 3, 5, 7]]
        first = build_system(CUBE, groups, [.75, .75])
        rename = lambda v: 17 + 3 * v
        faces = [[rename(v) for v in reversed(face)] for face in reversed(CUBE)]
        second = build_system(faces, [[rename(v) for v in g] for g in groups], [.75, .75])
        self.assertEqual([c["key"] for c in first["corner_types"]],
                         [c["key"] for c in second["corner_types"]])
        self.assertEqual({r["counts"] for r in first["fan_rows"]},
                         {r["counts"] for r in second["fan_rows"]})
        self.assertEqual(set(first["consequence_equations"]), set(second["consequence_equations"]))

    def test_rejects_open_complex_and_disconnected_vertex_link(self):
        with self.assertRaisesRegex(ValueError, "exactly two"):
            build_system([[0, 1, 2]], [[0, 1, 2]], [.375])
        faces = TETRA + [tuple(0 if v == 0 else v + 3 for v in f) for f in TETRA]
        with self.assertRaisesRegex(ValueError, "disconnected link"):
            build_system(faces, [list(range(7))], [.375])


class J88EncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = build_system(J88_FACES, J88_GROUPS, J88_RADII)

    def test_d_equation_and_c0_complement(self):
        model = self.model
        a, b, c, d = model["radius_variables"]
        expected_d = d ** 2 - (a + c) * d + (a - 1) * (c - 1)
        d_corner = next(x for x in model["corner_types"] if x["key"] == (3, 0, 2, 3))
        self.assertEqual(sp.expand(d_corner["N"] + expected_d), 0)
        single = next(r for r in model["sparse_relations"] if r["method"] == "single angle fixed exactly")
        self.assertFalse(sp.cancel(single["equations"][0] / expected_d).free_symbols)
        square = next(x for x in model["corner_types"] if x["key"] == (2, 0, 2, 4))
        triangle = next(x for x in model["corner_types"] if x["key"] == (2, 0, 3, 3))
        expected_c0 = square["N"] ** 2 * triangle["M"] - triangle["N"] ** 2 * square["M"]
        complement = next(r for r in model["sparse_relations"]
                          if r["rhs_pi"] == 1 and sum(x != 0 for x in r["counts"]) == 2)
        self.assertFalse(sp.cancel(expected_c0 / complement["equations"][0]).free_symbols)
        e = sp.Symbol("e")
        aa = (e ** 4 + 3 * e ** 2 - 1) / (4 * e ** 2 - 1)
        cc = aa + (1 - 2 * e ** 2) * (1 + e) / (4 * e ** 2 - 1)
        exact_family = {a: aa, c: cc, d: cc + e}
        self.assertEqual(sp.cancel(expected_d.subs(exact_family)), 0)
        self.assertEqual(sp.cancel(expected_c0.subs(exact_family)), 0)

    def test_distinct_fans_on_one_sphere_remain_distinct(self):
        fans = self.model["fan_counts_by_vertex"]
        self.assertEqual(fans[0], fans[5])
        self.assertEqual(fans[7], fans[10])
        self.assertNotEqual(fans[0], fans[7])
        self.assertEqual(len(self.model["fan_rows"]), 5)
        self.assertTrue(any(r["origin"] == "difference of vertex fans" and r["equations"]
                            for r in self.model["sparse_relations"]))

    def test_numerical_witness_obeys_full_original_equations(self):
        residual = max(abs(float(e.evalf(subs=self.model["witness"])))
                       for e in self.model["equations"])
        self.assertLess(residual, 1e-10)


class AngleIdentityTests(unittest.TestCase):
    @staticmethod
    def corner(angle):
        n = sp.cos(angle)
        return dict(N=n, M=sp.S.One, S=sp.sin(angle) ** 2,
                    D_left=sp.S.One, D_right=sp.S.One)

    def test_three_angle_consequences_with_signed_counts_and_odd_turns(self):
        cases = [((2, 2, 1), 1, (sp.pi / 6, sp.pi / 4, sp.pi / 6)),
                 ((2, -2, 1), 0, (sp.pi / 6, sp.pi / 3, sp.pi / 3)),
                 ((2, 2, -1), 0, (sp.pi / 8, sp.pi / 8, sp.pi / 2)),
                 ((2, -2, -1), -1, (sp.pi / 8, 3 * sp.pi / 8, sp.pi / 2))]
        for counts, rhs, angles in cases:
            with self.subTest(counts=counts, rhs=rhs):
                self.assertEqual(sp.simplify(sum(c * a for c, a in zip(counts, angles)) - rhs * sp.pi), 0)
                equations, _, _ = _relation_consequences(counts, sp.Integer(rhs), [self.corner(a) for a in angles])
                self.assertEqual(len(equations), 1)
                self.assertEqual(sp.simplify(sp.expand_trig(equations[0])), 0)
                self.assertFalse(equations[0].has(sp.Float))

    def test_two_angle_consequences_keep_exact_sign_when_possible(self):
        cases = [((1, 1), 1, (sp.pi / 4, 3 * sp.pi / 4)),
                 ((1, -1), 0, (sp.pi / 5, sp.pi / 5)),
                 ((4, 1), 2, (2 * sp.pi / 5, 2 * sp.pi / 5)),
                 ((1, -2), -1, (sp.pi / 3, 2 * sp.pi / 3))]
        for counts, rhs, angles in cases:
            equations, method, _ = _relation_consequences(counts, sp.Integer(rhs), [self.corner(a) for a in angles])
            self.assertIn("without squaring", method)
            self.assertEqual(sp.simplify(equations[0]), 0)
            self.assertFalse(equations[0].has(sp.Float))


if __name__ == "__main__":
    unittest.main()
