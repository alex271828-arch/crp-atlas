"""Local algebraic controls for the v2 reduction stages."""
import unittest
from unittest.mock import patch

import sympy as sp

from reduce import (_groebner_solutions, _quadratic_pair_consequences,
                    _radical_square_classes, solve_model)


class ReductionV2Tests(unittest.TestCase):
    def test_underconstrained_system_is_deferred_before_groebner(self):
        x, y = sp.symbols("x y")
        trace = []
        with patch.object(sp, "groebner", side_effect=AssertionError("must not run")):
            result = _groebner_solutions([x*y-1], [x, y], {x: 1, y: 1}, trace)
        self.assertEqual(result, [])
        self.assertEqual(trace[-1]["reason"], "structurally_underconstrained")

    def test_matching_detects_local_variable_deficiency(self):
        x, y, z = sp.symbols("x y z")
        trace = []
        with patch.object(sp, "groebner", side_effect=AssertionError("must not run")):
            result = _groebner_solutions([x+y-2, z-1, z*z-1], [x, y, z],
                                         {x: 1, y: 1, z: 1}, trace)
        self.assertEqual(result, [])
        self.assertEqual(trace[-1]["matching"], 2)

    def test_quadratic_cancellation_is_an_exact_ideal_combination(self):
        u, v = sp.symbols("u v")
        first, second = u*u-v, 2*u*u+u-3
        trace = []
        consequences = _quadratic_pair_consequences([first, second], [u, v], trace)
        self.assertIn(-u-2*v+3, consequences)
        self.assertEqual(sp.expand(consequences[0]-(2*first-second)), 0)
        self.assertEqual(sp.degree(consequences[0], u), 1)

    def test_radicals_share_only_exact_square_classes(self):
        x = sp.Symbol("x", positive=True)
        y0, y1 = sp.symbols("y0 y1", positive=True)
        radicands = [(x+1)**2*(x+2), 4*(x+2)/(x+3)**2]
        model = {"corner_types": [{"y": y0, "S": radicands[0]},
                                  {"y": y1, "S": radicands[1]}]}
        replacements, equations, auxiliaries = _radical_square_classes(
            model, {}, {x: 1.0, y0: 12**0.5, y1: 0.75**0.5}, [])
        self.assertEqual(len(auxiliaries), 1)
        q = auxiliaries[0]
        for y, radicand in zip((y0, y1), radicands):
            numerator = sp.cancel(replacements[y]**2-radicand).as_numer_denom()[0]
            self.assertEqual(sp.rem(numerator, equations[0], q), 0)
            self.assertFalse(replacements[y].has(sp.Float))

    def test_patch_equations_can_solve_before_full_fan(self):
        a, b = sp.symbols("s0 s1", positive=True)
        model = dict(radius_variables=(a, b), variables=(a, b),
                     patch_equations=[a-1, b*b-2],
                     witness={a: 1.0, b: 2**0.5}, equations=[])
        result = solve_model(model, max_seconds=2)
        self.assertEqual(result["status"], "ALGEBRAIC_CANDIDATE", result.get("reason"))
        self.assertEqual(len(result["phases"]), 1)
        self.assertEqual(result["phases"][0]["name"], "patch_and_radius_relations")

    def test_compact_positive_radical_closes_a_remaining_radius(self):
        a, b, y = sp.symbols("s0 s1 y0", positive=True)
        model = dict(radius_variables=(a, b), variables=(a, b, y),
                     patch_equations=[a-1], short_radical_equations=[y-2],
                     corner_types=[{"y": y, "S": b}],
                     witness={a: 1.0, b: 4.0, y: 2.0}, equations=[])
        result = solve_model(model, max_seconds=2)
        self.assertEqual(result["status"], "ALGEBRAIC_CANDIDATE", result.get("reason"))
        self.assertEqual(result["candidate"]["squared_radii"], [1, 4])
        self.assertEqual(result["phases"][-1]["name"], "compact_signed_radical_relations")


if __name__ == "__main__":
    unittest.main()
