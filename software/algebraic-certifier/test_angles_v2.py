"""Exact tests of short angle consequences and their retained sign branches."""
import unittest

import sympy as sp

try:
    from .encoder import (build_system, _relation_consequences,
                          _short_radical_relation, _bounded_short_relations)
except ImportError:
    from encoder import (build_system, _relation_consequences,
                         _short_radical_relation, _bounded_short_relations)


J87_FACES = [[6, 5, 10], [0, 2, 3], [5, 1, 7], [8, 6, 10], [9, 3, 8],
             [0, 4, 1], [0, 3, 4], [9, 8, 10], [3, 2, 8], [7, 1, 4],
             [7, 4, 9], [10, 7, 9], [8, 2, 6], [5, 7, 10], [0, 1, 2],
             [2, 1, 5, 6], [4, 3, 9]]
J87_GROUPS = [[1, 2], [0, 3, 4, 5, 6], [7, 8], [9, 10]]
J87_RADII = [.3336959211742356673, 1.0314576181363191673,
             1.0692829503133840764, 1.2985615632098495685]


def angle_corners(angles, norms=None):
    if norms is None:
        norms = [1] * len(angles)
    corners, assignment = [], {}
    for i, (angle, norm) in enumerate(zip(angles, norms)):
        norm = sp.sympify(norm)
        y = sp.Symbol(f"y{i}", positive=True)
        corners.append(dict(N=sp.sqrt(norm) * sp.cos(angle), M=norm,
                            S=norm * sp.sin(angle) ** 2, y=y,
                            D_left=norm, D_right=sp.S.One))
        assignment[y] = sp.sqrt(norm) * sp.sin(angle)
    return corners, assignment


class ShortAngleIdentityTests(unittest.TestCase):
    def test_three_angle_gram_for_signed_and_multiplied_angles(self):
        cases = [((1, 1, -1), 0, (sp.pi / 6, sp.pi / 3, sp.pi / 2)),
                 ((1, 1, 1), 1, (sp.pi / 3,) * 3),
                 ((1, -1, 1), 1, (2 * sp.pi / 3, sp.pi / 3, 2 * sp.pi / 3)),
                 ((1, -1, -1), -1, (sp.pi / 6, 2 * sp.pi / 3, sp.pi / 2)),
                 ((2, 1, -1), 0, (sp.pi / 6, sp.pi / 6, sp.pi / 2))]
        for counts, rhs, angles in cases:
            for norms in [(1, 1, 1), (2, 3, 5)]:
                with self.subTest(counts=counts, rhs=rhs, norms=norms):
                    self.assertEqual(sp.simplify(sum(k * a for k, a in zip(counts, angles)) - rhs * sp.pi), 0)
                    corners, assignment = angle_corners(angles, norms)
                    equations, method, _ = _relation_consequences(counts, sp.Integer(rhs), corners)
                    self.assertIn("cosine Gram determinant", method)
                    self.assertEqual(len(equations), 1)
                    self.assertEqual(sp.simplify(equations[0]), 0)
                    mixed, signs, _ = _short_radical_relation(counts, sp.Integer(rhs), corners)
                    self.assertTrue(mixed)
                    self.assertTrue(all(sp.simplify(e.subs(assignment)) == 0 for e in mixed))
                    self.assertTrue(all(sp.simplify(e.subs(assignment)).is_positive for e in signs))
                    self.assertFalse(any(e.has(sp.Float) for e in equations + mixed + signs))

    def test_four_angle_pair_norm_and_signed_mixed_closure(self):
        cases = [((1, 1, 1, 1), 2, (sp.pi / 3, sp.pi / 2, sp.pi / 2, 2 * sp.pi / 3)),
                 ((1, 1, -1, -1), 0, (sp.pi / 6, sp.pi / 3, sp.pi / 4, sp.pi / 4)),
                 ((1, 1, 1, 1), 1, (sp.pi / 6, sp.pi / 6, sp.pi / 3, sp.pi / 3))]
        for counts, rhs, angles in cases:
            corners, assignment = angle_corners(angles, [2, 2, 3, 3])
            equations, method, _ = _relation_consequences(counts, sp.Integer(rhs), corners)
            self.assertIn("four-angle", method)
            self.assertEqual(sp.simplify(equations[0]), 0)
            mixed, signs, _ = _short_radical_relation(counts, sp.Integer(rhs), corners)
            self.assertTrue(all(sp.simplify(e.subs(assignment)) == 0 for e in mixed))
            self.assertTrue(all(sp.simplify(e.subs(assignment)).is_positive for e in signs))

    def test_wrong_radical_sign_is_detected_by_retained_equations(self):
        corners, assignment = angle_corners((sp.pi / 6, sp.pi / 3, sp.pi / 2))
        equations, _, _ = _short_radical_relation((1, 1, -1), sp.S.Zero, corners)
        self.assertTrue(all(sp.simplify(e.subs(assignment)) == 0 for e in equations))
        assignment[corners[0]["y"]] *= -1
        self.assertTrue(any(sp.simplify(e.subs(assignment)) != 0 for e in equations))

    def test_sparse_combinations_stay_in_exact_fan_row_space(self):
        rows = [dict(counts=row) for row in [(1, 1, 0, 1, 0, 0),
                                             (0, 1, 1, 0, 1, 0),
                                             (0, 0, 1, 1, 1, 1)]]
        found, diagnostics = _bounded_short_relations(rows)
        reverse, _ = _bounded_short_relations(list(reversed(rows)))
        self.assertEqual(found, reverse)
        self.assertIn(((1, 0, 0, 2, 0, 1), sp.Integer(2)), found)
        matrix = sp.Matrix([list(row["counts"]) + [2] for row in rows])
        for (counts, rhs) in found:
            self.assertEqual(matrix.col_join(sp.Matrix([list(counts) + [rhs]])).rank(), matrix.rank())
        self.assertFalse(diagnostics["budget_exhausted"])
        _, bounded = _bounded_short_relations(rows, max_combinations=2)
        self.assertEqual(bounded["attempted_combinations"], 2)
        self.assertTrue(bounded["budget_exhausted"])

    def test_large_radius_only_expansion_is_kept_as_mixed_equations(self):
        variables = sp.symbols("s0:4", positive=True)
        corners = [dict(N=s, M=sp.S.One, S=1 - s * s,
                        D_left=sp.S.One, D_right=sp.S.One,
                        y=sp.Symbol(f"y{i}", positive=True)) for i, s in enumerate(variables)]
        equations, _, _ = _relation_consequences((1, 1, 1, 1), sp.Integer(2), corners)
        self.assertEqual(equations, [])
        mixed, _, _ = _short_radical_relation((1, 1, 1, 1), sp.Integer(2), corners)
        self.assertEqual(len(mixed), 2)
        self.assertTrue(all(sp.Poly(e, *variables, *(c["y"] for c in corners)).total_degree() <= 2 for e in mixed))


class J87ShortAngleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = build_system(J87_FACES, J87_GROUPS, J87_RADII)

    def test_manual_patch_quadratic_is_generated_from_fan_row_space(self):
        model = self.model
        a, b, _, _ = model["radius_variables"]
        expected = 3 * (a - b) ** 2 - 4 * (a + b) + 4
        matching = [r for r in model["sparse_relations"]
                    if sum(abs(k) for k in r["counts"]) == 3
                    and r["rhs_pi"] == 0]
        self.assertEqual(len(matching), 1)
        relation = matching[0]
        self.assertIn("without squaring", relation["method"])
        self.assertEqual(sp.div(relation["equations"][0], expected, a, b)[1], 0)
        y5 = model["corner_types"][5]["y"]
        y8 = model["corner_types"][8]["y"]
        exact_mixed = 2 * b * (3 * b - a) - y5 * y8
        self.assertTrue(any(sp.expand(e - exact_mixed) == 0 for e in relation["short_radical_equations"]))

    def test_all_short_equations_match_the_original_positive_branch(self):
        equations = self.model["short_radical_equations"]
        self.assertGreaterEqual(len(equations), 10)
        self.assertFalse(any(e.has(sp.Float) for e in equations))
        residual = max(abs(float(e.evalf(subs=self.model["witness"]))) for e in equations)
        self.assertLess(residual, 1e-10)

    def test_relabeling_does_not_change_short_relations(self):
        rename = lambda vertex: 20 + 7 * vertex
        faces = [[rename(v) for v in reversed(face)] for face in reversed(J87_FACES)]
        groups = [[rename(v) for v in group] for group in J87_GROUPS]
        second = build_system(faces, groups, J87_RADII)
        first_keys = {(r["counts"], r["rhs_pi"]) for r in self.model["sparse_relations"]}
        second_keys = {(r["counts"], r["rhs_pi"]) for r in second["sparse_relations"]}
        self.assertEqual(first_keys, second_keys)
        self.assertEqual(set(self.model["short_radical_equations"]), set(second["short_radical_equations"]))


if __name__ == "__main__":
    unittest.main()
