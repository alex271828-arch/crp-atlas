import test_support as reference

"""Independent arithmetic checks and known-target pruning regressions."""
import json, math, random, unittest
from fractions import Fraction
from pathlib import Path
import numpy as np
import certified_geometry as C
import link_growth as G
import surface_pruning as P
import boundary_filter as B
import census_alphabet as E


class PruningTests(unittest.TestCase):
    def test_integer_rounding(self):
        rng = random.Random(111)
        for _ in range(1000):
            a, b = sorted(rng.sample(range(-2 * C.S, 2 * C.S), 2))
            c, d = sorted(rng.sample(range(-2 * C.S, 2 * C.S), 2))
            x, y = C.I(a, b), C.I(c, d)
            for op in ("mul", "div"):
                if op == "div" and c <= 0 <= d:
                    continue
                z = x * y if op == "mul" else x / y
                vals = [
                    Fraction(i * j, C.S**2) if op == "mul" else Fraction(i, j)
                    for i in (a, b)
                    for j in (c, d)
                ]
                self.assertLessEqual(Fraction(z.a, C.S), min(vals))
                self.assertGreaterEqual(Fraction(z.b, C.S), max(vals))
            z = C.sqrt(C.I(abs(a), max(abs(a), abs(b))))
            self.assertLessEqual(z.a * z.a, abs(a) * C.S)
            self.assertGreaterEqual(z.b * z.b, max(abs(a), abs(b)) * C.S)

    def test_engine_inserts_mandatory_face_before_branching(self):
        import link_engine as LE

        state = LE.initial(((0, 1, 2), (0, 2, 3), (0, 3, 1)), ((3, 3, 3),))
        found = []
        LE.run(state, 10, lambda fs, code: found.append(fs))
        self.assertTrue(state["exhausted"])
        self.assertEqual(len(found), 1)
        self.assertEqual(len(found[0]), 4)
        self.assertEqual(state["stats"]["faces_generated"], 0)

    def test_cached_transfers_match_raw(self):
        rng = random.Random(114)
        fours = [w for w in E.SIGS if len(w) == 4]
        fives = [w for w in E.SIGS if len(w) == 5]

        def outcome(fn, *args):
            try:
                r = fn(*args)
                return (
                    {k: (v.a, v.b) for k, v in r.items()}
                    if isinstance(r, dict)
                    else tuple((v.a, v.b) for v in r)
                )
            except C.Impossible:
                return None

        for _ in range(500):
            w = rng.choice(fours)
            k = rng.randrange(4)
            a, b = sorted(rng.sample(range(-C.S, C.S + 1), 2))
            x = C.I(a, b)
            self.assertEqual(outcome(C.quad, w, k, x), outcome(C._quad_raw, w, k, x))
        for _ in range(150):
            w = rng.choice(fives)
            a, b = sorted(rng.sample(range(-C.S, C.S + 1), 2))
            c, d = sorted(rng.sample(range(-C.S, C.S + 1), 2))
            x, y = C.I(a, b), C.I(c, d)
            self.assertEqual(outcome(C.pent, w, x, y), outcome(C._pent_raw, w, x, y))

    def test_equal_side_singular_limit(self):
        for n in (3, 4, 5, 6, 8, 10):
            z = C.tri(C.SIDE[n], C.SIDE[n], C.I(C.S))
            self.assertLessEqual(z.a, 0)
            self.assertGreaterEqual(z.b, 0)
            self.assertLess(z.width(), 100)

    def test_monotone_angle_bounds(self):
        import mpmath as mp

        rng = random.Random(113)
        with mp.workdps(60):
            for _ in range(250):
                a, b = sorted(rng.sample(range(-C.S, C.S + 1), 2))
                c, d = sorted(rng.sample(range(-C.S, C.S + 1), 2))
                x, y = C.I(a, b), C.I(c, d)
                for fn, sign in ((C.angle_sum, 1), (C.angle_difference, -1)):
                    vals = []
                    for u in (a, (a + b) // 2, b):
                        for v in (c, (c + d) // 2, d):
                            angle = mp.acos(mp.mpf(u) / C.S) + sign * mp.acos(
                                mp.mpf(v) / C.S
                            )
                            if 0 <= angle <= mp.pi:
                                vals.append(mp.cos(angle))
                    try:
                        z = fn(x, y)
                    except C.Impossible:
                        self.assertFalse(vals)
                        continue
                    for val in vals:
                        self.assertLessEqual(mp.mpf(z.a) / C.S, val)
                        self.assertGreaterEqual(mp.mpf(z.b) / C.S, val)

    def test_run_dp_against_reference(self):
        rng = random.Random(112)
        pool = [(a, b, 0, (a, 3, b)) for a in B.ANGLE for b in B.ANGLE] + [
            (a, a, 0, (a,)) for a in B.ANGLE
        ]
        for _ in range(1000):
            n = rng.randrange(3, 12)
            ds = {v: tuple(rng.sample(pool, rng.randrange(1, 8))) for v in range(n)}
            rows = tuple(
                tuple(sorted({(a, b, len(g)) for a, b, k, g in ds[v]}))
                for v in range(n)
            )
            self.assertEqual(
                P.run_possible(rows),
                reference.cycle_minimum(list(range(n)), ds, True) < 10**9,
            )

    def test_singular_triangle_undivided_equation(self):
        with self.assertRaises(C.Impossible):
            C.tri(C.I(0), C.I(C.S), C.I(C.S))
        x = C.tri(C.I(C.S), C.I(C.S), C.I(C.S))
        self.assertEqual((x.a, x.b), (-C.S, C.S))

    def test_all_catalog_geometry(self):
        count = 0
        for s in json.loads(Path("testdata/crp_catalog.json").read_text())["solids"]:
            if any(len(f) not in C.SIDE for f in s["faces"]):
                continue
            fs = s["faces"]
            x = np.array(s["vertices"])
            norm = []
            ef = {}
            for i, f in enumerate(fs):
                n = np.cross(x[f[1]] - x[f[0]], x[f[2]] - x[f[1]])
                n /= np.linalg.norm(n)
                norm.append(n)
                for j, v in enumerate(f):
                    ef.setdefault(frozenset((v, f[j - 1])), []).append(i)
            box = {
                e: C.cosine(
                    C.I(
                        math.floor((-np.dot(norm[a], norm[b]) - 1e-6) * C.S),
                        math.ceil((-np.dot(norm[a], norm[b]) + 1e-6) * C.S),
                    )
                )
                for e, (a, b) in ef.items()
            }
            info = G.analyze(fs, tuple(B.WORDS))
            vw = C.vertex_words(info)
            order, _ = C.plan(vw)
            with self.subTest(solid=s["primary_name"]):
                C.evaluate(vw, order, box)
                C.convex_check(fs, box)
            count += 1
        self.assertEqual(count, 120)

    def test_j40_pinched_boundary(self):
        d = json.loads(Path("testdata/pass20_pivot_dead_end.json").read_text())
        fs = frozenset(map(tuple, d["faces"]))
        new = tuple(d["new"][0])
        menu = tuple(sorted(E.SIGS[E.BASELINE_PASS[20] :]))
        info = G.analyze(fs, menu)
        self.assertFalse(P.propagate(fs, info, menu)["excluded"])
        target = G.canonical_code(fs | {new}, d["pivot"])
        self.assertIn(
            target,
            {
                G.canonical_code(fs | {f}, d["pivot"])
                for f in G.next_faces(fs, d["pivot"], menu, info)
            },
        )
        info = G.analyze(fs | {new}, menu)
        self.assertTrue(info["pinched"])
        self.assertFalse(P.propagate(fs | {new}, info, menu)["excluded"])
        self.assertFalse(C.solve(fs | {new}, info, time_cap=0.02)[0])

    def test_timeout_preserves_box(self):
        fs = ((0, 1, 2), (0, 2, 3), (0, 3, 1))
        info = G.analyze(fs, tuple(B.WORDS))
        excluded, region, count = C.solve(fs, info, time_cap=0)
        self.assertFalse(excluded)
        self.assertEqual(count, 0)
        self.assertTrue(region["boxes"])
        self.assertFalse(C.solve(fs, info, region, time_cap=0.1)[0])

    def test_unresolved_box_does_not_hide_other_boxes(self):
        from unittest.mock import patch

        fs = ((0, 1, 2), (0, 2, 3), (0, 3, 1))
        info = G.analyze(fs, tuple(B.WORDS))
        e = frozenset((0, 1))
        parent = {"faces": fs, "boxes": [{e: (C.S, C.S)}, {e: (0, 0)}]}

        def evaluate(vw, order, box):
            if box[e].a == C.S:
                raise C.Impossible
            return box

        with patch.object(C, "plan", return_value=([], [e])), patch.object(
            C, "evaluate", side_effect=evaluate
        ), patch.object(C, "convex_check"):
            excluded, region, count = C.solve(fs, info, parent, time_cap=1)
        self.assertFalse(excluded)
        self.assertEqual(count, 2)
        self.assertEqual(region["boxes"], [{e: (0, 0)}])

    def test_known_target_growth_with_filters(self):
        steps = pinches = 0
        rng = random.Random(111)
        for s in json.loads(Path("testdata/crp_catalog.json").read_text())["solids"]:
            if any(len(f) not in C.SIDE for f in s["faces"]):
                continue
            full = frozenset(map(tuple, s["faces"]))
            menu = tuple(B.WORDS)
            import rim_geometry as R

            coords = np.array(s["vertices"])
            norm = []
            ef = {}
            for fi, f in enumerate(s["faces"]):
                nv = np.cross(coords[f[1]] - coords[f[0]], coords[f[2]] - coords[f[1]])
                nv /= np.linalg.norm(nv)
                norm.append(nv)
                for j, v in enumerate(f):
                    ef.setdefault(frozenset((v, f[j - 1])), []).append(fi)
            angles = {
                e: C.cosine(
                    C.I(
                        math.floor((-np.dot(norm[a], norm[b]) - 1e-6) * C.S),
                        math.ceil((-np.dot(norm[a], norm[b]) + 1e-6) * C.S),
                    )
                )
                for e, (a, b) in ef.items()
            }
            # Random edge-connected face prefixes exercise disconnected fans.
            for trial in range(2):
                fs = frozenset({rng.choice(sorted(full))})
                while True:
                    info = G.analyze(fs, menu)
                    self.assertIsNotNone(info)
                    with self.subTest(solid=s["primary_name"], faces=len(fs)):
                        verdict = P.propagate(fs, info, menu)
                        self.assertFalse(verdict["excluded"])
                        canon = lambda f: min(
                            tuple(f[i:]) + tuple(f[:i]) for i in range(len(f))
                        )
                        self.assertLessEqual(
                            {canon(f) for f in verdict["faces"]},
                            {canon(f) for f in full},
                        )
                        R.check(R.fans(info), angles)
                    pinches += bool(info["pinched"])
                    steps += 1
                    if fs == full:
                        break
                    choices = [
                        f
                        for f in full - fs
                        if any(
                            (f[i], f[i - 1]) in info["directed"] for i in range(len(f))
                        )
                    ]
                    fs = fs | {rng.choice(sorted(choices))}
        self.assertGreater(steps, 2000)
        self.assertGreater(pinches, 20)


if __name__ == "__main__":
    unittest.main()
