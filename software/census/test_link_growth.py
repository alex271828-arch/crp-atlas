import itertools
import json
import random
import tempfile
import unittest
import test_support as reference
from pathlib import Path
import boundary_filter as BF
import census_alphabet as E
import link_growth as G
import link_engine as LE


def rotations(word):
    return {t[i:] + t[:i] for t in (word, word[::-1]) for i in range(len(word))}


def reference_support(paths, menu):
    result = []
    for sig in menu:
        missing = len(sig) - sum(map(len, paths))
        good = False
        for order in itertools.permutations(paths):
            for gaps in itertools.product(
                range(1, max(1, missing) + 1), repeat=len(paths)
            ):
                if sum(gaps) != missing:
                    continue
                pattern = tuple(x for p, n in zip(order, gaps) for x in p + (None,) * n)
                if any(
                    all(a is None or a == b for a, b in zip(pattern, w))
                    for w in rotations(sig)
                ):
                    good = True
                    break
            if good:
                break
        if good:
            result.append(sig)
    return tuple(result)


class LinkTests(unittest.TestCase):
    def test_joint_domains_against_independent_gap_enumerator(self):
        menu = tuple(BF.WORDS)
        paths = [(n,) for n in G.ANGLE] + list(itertools.product(G.ANGLE, repeat=2))
        for a in paths:
            self.assertEqual(G.supported((a,), menu), reference_support((a,), menu))
            for b in paths:
                ps = tuple(sorted((a, b)))
                self.assertEqual(G.supported(ps, menu), reference_support(ps, menu), ps)

    def test_link_structure(self):
        menu = tuple(BF.WORDS)
        self.assertIsNotNone(G.vertex_domain(((1, 2, 3), (3, 4, 3)), menu))
        self.assertIsNone(G.vertex_domain(((1, 2, 3), (1, 3, 3)), menu))
        self.assertIsNone(
            G.vertex_domain(((1, 2, 3), (2, 3, 3), (3, 1, 3), (4, 5, 3)), menu)
        )
        self.assertIsNone(G.vertex_domain(((1, 2, 6), (2, 3, 6), (3, 1, 6)), menu))
        # Two disjoint two-face fans cannot fit a degree-five link: two gaps needed.
        self.assertIsNone(
            G.vertex_domain(((1, 2, 3), (2, 3, 3), (4, 5, 3), (5, 6, 3)), menu)
        )

    def test_face_enumeration_against_brute_force(self):
        menu = tuple(w for w in BF.WORDS if set(w) <= {3, 4})
        for faces in (
            frozenset({(0, 1, 2)}),
            frozenset({(0, 1, 2, 3)}),
            frozenset({(0, 1, 2), (0, 3, 4), (2, 1, 4, 3)}),
        ):
            state = G.analyze(faces, menu)
            pivot = 0
            u = min(p[1] for p in state["domains"][pivot][1])
            base = max(state["vertices"]) + 1
            expected = set()
            for n in (3, 4):
                pool = sorted(state["opened"]) + list(range(base, base + n - 2))
                for tail in itertools.permutations(
                    [x for x in pool if x not in (u, pivot)], n - 2
                ):
                    fresh = [x for x in tail if x >= base]
                    if fresh != list(range(base, base + len(fresh))):
                        continue
                    f = (u, pivot) + tail
                    if BF.intersection_obstruction(list(faces) + [f]):
                        continue
                    if G.analyze(faces | {f}, menu) is not None:
                        expected.add(f)
            actual = {
                f
                for f in G.next_faces(faces, pivot, menu, state)
                if G.analyze(faces | {f}, menu) is not None
            }
            self.assertEqual(actual, expected)

    def test_j40_pinch_is_retained_without_switching_pivot(self):
        d = json.loads(
            Path(__file__)
            .with_name("testdata")
            .joinpath("pass20_pivot_dead_end.json")
            .read_text()
        )
        b = frozenset(map(tuple, d["faces"]))
        v = d["pivot"]
        menu = tuple(sorted(E.SIGS[E.BASELINE_PASS[20] :]))
        face = tuple(d["new"][0])
        target = G.canonical_code(b | {face}, v)
        self.assertIn(
            target, {G.canonical_code(b | {f}, v) for f in G.next_faces(b, v, menu)}
        )
        info = G.analyze(b | {face}, menu)
        self.assertNotIn(v, info["opened"])
        self.assertIn(32, info["pinched"])

    def test_canonicalization_and_root(self):
        fs = ((0, 1, 2), (2, 1, 3, 4))
        code = G.canonical_code(fs, 0)
        rng = random.Random(110)
        for _ in range(12):
            labels = list(range(20, 25))
            rng.shuffle(labels)
            rel = dict(enumerate(labels))
            other = [tuple(rel[v] for v in f[::-1]) for f in reversed(fs)]
            self.assertEqual(code, G.canonical_code(other, rel[0]))
        self.assertNotEqual(code, G.canonical_code(fs, 1))

    def test_known_host_next_faces(self):
        cat = json.loads(
            Path(__file__)
            .with_name("testdata")
            .joinpath("crp_catalog.json")
            .read_text()
        )["solids"]
        tested = 0
        for solid in cat:
            if not solid["primary_name"].startswith(("J40 ", "J41 ")):
                continue
            full = frozenset(map(tuple, solid["faces"]))
            menu = tuple(sorted(E.SIGS[E.BASELINE_PASS[20] :]))
            start = next(
                v
                for v in reference.vertices_of(full)
                if reference.link_signature(full, v) == E.SIGS[E.BASELINE_PASS[20]]
            )
            b = frozenset(f for f in full if start in f)
            while b != full:
                info = G.analyze(b, menu)
                self.assertIsNotNone(info)
                pivot = G.choose_pivot(info)
                while pivot in info["opened"]:
                    u = min(p[1] for p in info["domains"][pivot][1])
                    target = next(
                        f
                        for f in full - b
                        if any(
                            f[i] == u and f[(i + 1) % len(f)] == pivot
                            for i in range(len(f))
                        )
                    )
                    child = b | {target}
                    code = G.canonical_code(child, pivot)
                    self.assertIn(
                        code,
                        {
                            G.canonical_code(b | {f}, pivot)
                            for f in G.next_faces(b, pivot, menu, info)
                        },
                    )
                    b = child
                    info = G.analyze(b, menu)
                    self.assertIsNotNone(info)
                    tested += 1
        self.assertGreater(tested, 50)


class EngineTests(unittest.TestCase):
    def enumerate(self, seed, menu, split=False):
        state = LE.initial(seed, menu)
        results = set()
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "state.pkl"
            for _ in range(10000):
                LE.run(
                    state,
                    10,
                    lambda f, c: results.add(c),
                    path,
                    max_steps=2 if split else 100000,
                )
                if state["exhausted"]:
                    break
                state = LE.load(path, menu)
            self.assertTrue(state["exhausted"])
        return results, state

    def test_tetrahedron_cube_octahedron(self):
        cat = json.loads(
            Path(__file__)
            .with_name("testdata")
            .joinpath("crp_catalog.json")
            .read_text()
        )["solids"]
        for word in ((3, 3, 3), (4, 4, 4), (3, 3, 3, 3)):
            solid = next(
                s
                for s in cat
                if all(
                    reference.link_signature(frozenset(map(tuple, s["faces"])), v)
                    == word
                    for v in reference.vertices_of(s["faces"])
                )
            )
            fs = frozenset(map(tuple, solid["faces"]))
            v = min(reference.vertices_of(fs))
            seed = [f for f in fs if v in f]
            found, state = self.enumerate(seed, (word,))
            self.assertEqual(found, {G.canonical_code(fs)})

    def test_resume_in_middle_of_link(self):
        seed = ((0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 1))
        menu = ((3, 3, 3, 3),)
        probe = LE.initial(seed, menu)
        LE.run(probe, 10, lambda f, c: None, max_steps=1)
        self.assertTrue(any(a[0] == "inner" for a in probe["agenda"]))
        a, sa = self.enumerate(seed, menu)
        b, sb = self.enumerate(seed, menu, True)
        self.assertEqual(len(a), 1)
        self.assertEqual(a, b)
        self.assertEqual(sa["seen_inner"], sb["seen_inner"])
        self.assertEqual(sa["seen_outer"], sb["seen_outer"])

    def test_closed_candidate_is_retained_without_numerical_gate(self):
        seed = ((0, 1, 2), (0, 2, 3), (0, 3, 1))
        menu = ((3, 3, 3),)
        found, _ = self.enumerate(seed, menu)
        self.assertEqual(len(found), 1)

    def test_checkpoint_menu_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "s.pkl"
            LE.save(LE.initial(((0, 1, 2),), ((3, 3, 3),)), p)
            with self.assertRaises(ValueError):
                LE.load(p, ((4, 4, 4),))


if __name__ == "__main__":
    unittest.main()
