import test_support as reference

"""New-alphabet coverage, exact angles, and large-face end-to-end regressions."""
import contextlib, io, itertools, json, math, tempfile, unittest
from fractions import Fraction
from pathlib import Path
import numpy as np
import mpmath as mp
import census_alphabet as A
import face_types as F
import boundary_filter as B
import certified_geometry as C
import link_growth as G
import link_engine as LE
import rim_geometry as R

HERE = Path(__file__).resolve().parent
CAT = json.loads((HERE / "testdata/crp_catalog.json").read_text())["solids"]


def prism(n):
    return next(
        s
        for s in CAT
        if len(s["vertices"]) == 2 * n
        and sorted(map(len, s["faces"])) == sorted([n, n] + [4] * n)
    )


def coplanar_large_face():
    d = json.loads((HERE / "testdata" / "large_face_coplanar.json").read_text())
    return [tuple(f) for f in d["faces"]], np.array(d["vertices"])


def known_angles(fs, x):
    normals = []
    ef = {}
    for fi, f in enumerate(fs):
        v = np.cross(x[f[1]] - x[f[0]], x[f[2]] - x[f[0]])
        v /= np.linalg.norm(v)
        normals.append(v)
        for i, a in enumerate(f):
            ef.setdefault(frozenset((a, f[i - 1])), []).append(fi)
    return {
        e: C.cosine(
            C.I(
                math.floor((-float(np.dot(normals[i], normals[j])) - 1e-6) * C.S),
                math.ceil((-float(np.dot(normals[i], normals[j])) + 1e-6) * C.S),
            )
        )
        for e, (i, j) in ef.items()
    }


class AlphabetTests(unittest.TestCase):
    def test_independent_rational_link_inventory(self):
        words = {}
        for d in range(3, 6):
            for w in itertools.product(F.FACE_SIZES, repeat=d):
                if sum((Fraction(n - 2, n) for n in w), Fraction()) < 2:
                    words[B.canonical(w)] = w
        self.assertEqual(set(words), set(A.SIGS))
        self.assertEqual(set(B.WORDS), set(A.SIGS))
        self.assertEqual(
            (len(A.NEW_LINKS), len(A.OLD_LINKS), len(A.SIGS)), (24, 39, 63)
        )
        self.assertEqual(F.ANGLE[8], 45)
        self.assertEqual(F.ANGLE[10], 48)
        self.assertNotIn((4, 8, 8), A.SIGS)
        self.assertEqual(
            min(F.FULL_TURN - sum(F.ANGLE[n] for n in w) for w in A.SIGS), 2
        )

    def test_old_pass_order_and_partition(self):
        old = json.loads((HERE / "provenance" / "v1_11_pass_order.json").read_text())
        self.assertEqual(A.ORDER[A.BASELINE_OFFSET :], old)
        for s in CAT:
            fs = frozenset(map(tuple, s["faces"]))
            types = {reference.link_signature(fs, v) for v in reference.vertices_of(fs)}
            k = min(A.SIGS.index(w) for w in types)
            self.assertEqual(k < 24, any(len(f) > 6 for f in fs))
            self.assertLessEqual(types, set(A.SIGS[k:]))

    def test_all_seed_stars(self):
        for k, w in enumerate(A.SIGS):
            fs = A.SEEDS[A.ORDER[k]]["central"]
            info = G.analyze(fs, tuple(A.SIGS[k:]))
            self.assertIsNotNone(info)
            self.assertFalse(info["domains"][0][0] is False)
            self.assertEqual(reference.link_signature(frozenset(map(tuple, fs)), 0), w)

    def test_exact_face_and_rim_constants(self):
        with mp.workdps(90):

            def enclosed(i, x):
                self.assertLessEqual(mp.mpf(i.a) / C.S, x + mp.mpf("1e-80"))
                self.assertGreaterEqual(mp.mpf(i.b) / C.S, x - mp.mpf("1e-80"))

            for n in F.FACE_SIZES:
                enclosed(C.SIDE[n], mp.cos(mp.pi * (n - 2) / n))
                enclosed(C.TURN[n], mp.cos(2 * mp.pi / n))
            for k, c in enumerate(R.COS):
                enclosed(c, mp.cos(mp.pi * k / 60))

    def test_large_face_generation_and_completion(self):
        for n in (8, 10):
            s = prism(n)
            full = frozenset(map(tuple, s["faces"]))
            missing = next(f for f in full if len(f) == n)
            fs = full - {missing}
            menu = ((4, 4, n),)
            info = G.analyze(fs, menu)
            pivot = min(info["opened"])
            target = G.canonical_code(full, pivot)
            found = any(
                G.canonical_code(fs | {f}, pivot) == target
                for f in G.next_faces(fs, pivot, menu, info)
            )
            self.assertTrue(found)
            state = LE.initial(fs, menu)
            closed = []
            LE.run(state, 10, lambda f, c: closed.append(c))
            self.assertTrue(state["exhausted"])
            self.assertEqual(closed, [G.canonical_code(full)])

    def test_coplanar_decagonal_body_is_retained(self):
        fs, x = coplanar_large_face()
        info = G.analyze(fs, tuple(A.SIGS))
        self.assertIsNotNone(info)
        a = known_angles(fs, x)
        vw = C.vertex_words(info)
        order, _ = C.plan(vw)
        C.evaluate(vw, order, a)
        C.convex_check(fs, a)
        R.check(R.fans(info), a)


if __name__ == "__main__":
    unittest.main()
