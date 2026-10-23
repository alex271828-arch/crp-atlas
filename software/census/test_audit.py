"""Independent combinatorial oracles and crash/restart regressions for v1.16."""

import contextlib
import io
import itertools
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import boundary_filter as B
import census
import certified_geometry as C
import link_engine as E
import link_growth as G
import surface_pruning as P
from storage import append, lock, rows

HERE = Path(__file__).resolve().parent
SEED = ((0, 1, 2), (0, 2, 3), (0, 3, 1))
MENU = ((3, 3, 3),)


class CombinatorialAudit(unittest.TestCase):
    def test_cyclic_support_against_exhaustive_assignments(self):
        rng = random.Random(116)
        pool = [(a, b, 0, (a, b)) for a in (3, 4, 5) for b in (3, 4, 5)]
        for _ in range(200):
            size = rng.randrange(2, 7)
            domains = {
                i: tuple(rng.sample(pool, rng.randrange(1, 4))) for i in range(size)
            }
            expected = {i: set() for i in domains}
            for assignment in itertools.product(*domains.values()):
                if all(
                    assignment[i][1] == assignment[(i + 1) % size][0] for i in domains
                ):
                    for i, row in enumerate(assignment):
                        expected[i].add(row)
            actual = B.supported_domains([list(domains)], domains)
            self.assertEqual({i: set(r) for i, r in actual.items()}, expected)

    def test_polygon_runs_against_exhaustive_assignments(self):
        rng = random.Random(117)
        pool = [(n, n, 1) for n in (3, 4, 5)] + [
            (a, b, 2) for a in (3, 4, 5) for b in (3, 4, 5)
        ]
        for _ in range(200):
            size = rng.randrange(3, 8)
            domains = tuple(
                tuple(rng.sample(pool, rng.randrange(1, 4))) for _ in range(size)
            )
            possible = False
            for assignment in itertools.product(*domains):
                if any(
                    assignment[i][1] != assignment[(i + 1) % size][0]
                    for i in range(size)
                ):
                    continue
                breaks = [i for i, row in enumerate(assignment) if row[2] > 1]
                if not breaks:
                    possible = assignment[0][0] == size
                else:
                    # A run begins at a break's outgoing edge and ends at the
                    # next break's incoming edge; it uses fewer than m edges.
                    possible = all(
                        ((b - a) % size or size) < assignment[a][1]
                        for a, b in zip(breaks, breaks[1:] + breaks[:1])
                    )
                if possible:
                    break
            self.assertEqual(P.run_possible(domains), possible)

    def test_canonical_code_invariant_under_face_rotations(self):
        rng = random.Random(118)
        catalog = json.loads((HERE / "testdata/crp_catalog.json").read_text())["solids"]
        for solid in catalog:
            faces = [tuple(f) for f in solid["faces"]]
            vertices = sorted({v for f in faces for v in f})
            shuffled = vertices[:]
            rng.shuffle(shuffled)
            relabel = dict(zip(vertices, shuffled))
            changed = []
            for f in faces:
                k = rng.randrange(len(f))
                changed.append(tuple(relabel[v] for v in (f[k:] + f[:k])[::-1]))
            rng.shuffle(changed)
            for root in (None, vertices[0]):
                self.assertEqual(
                    G.canonical_code(faces, root),
                    G.canonical_code(changed, None if root is None else relabel[root]),
                )

    def test_face_intersection_checked_incrementally(self):
        faces = [(0, 1, 2), (2, 1, 3, 4)]
        for new in [(0, 2, 4), (0, 1, 3), (2, 4, 5), (4, 3, 5)]:
            self.assertEqual(
                bool(B.intersection_obstruction(faces + [new])),
                any(B.intersection_obstruction([old, new]) for old in faces),
            )


class RecoveryAudit(unittest.TestCase):
    def test_journal_torn_tail_replayed(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "journal.jsonl"
            append(path, {"number": 1})
            with path.open("ab") as handle:
                handle.write(b'{"number": 2')
            with self.assertRaises(ValueError):
                rows(path)
            self.assertEqual(rows(path, repair_tail=True), [{"number": 1}])
            append(path, {"number": 2})
            self.assertEqual(rows(path), [{"number": 1}, {"number": 2}])

    def test_malformed_complete_record_is_not_discarded(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "journal.jsonl"
            path.write_bytes(b"{broken}\n")
            with self.assertRaises(json.JSONDecodeError):
                rows(path, repair_tail=True)
            self.assertEqual(path.read_bytes(), b"{broken}\n")

    def test_sink_failure_preserves_unexpanded_item(self):
        state = E.initial(SEED, MENU)
        original = list(state["agenda"])

        def fail(*_):
            raise OSError("Simulated journal failure")

        with self.assertRaises(OSError):
            E.run(state, 10, fail)
        self.assertEqual(state["agenda"], original)
        self.assertFalse(state["seen_outer"])
        found = []
        E.run(state, 10, lambda f, c: found.append(c))
        self.assertTrue(state["exhausted"])
        self.assertEqual(len(found), 1)

    def test_partial_generation_failure_keeps_all_siblings_replayable(self):
        seed = ((0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 1))
        menu = ((3, 3, 3, 3),)
        state = E.initial(seed, menu)
        E.run(state, 10, lambda *_: None, max_steps=1)
        original = list(state["agenda"])
        original_seen = set(state["seen_inner"])
        generator = G.next_faces

        def interrupted(*args):
            yield next(generator(*args))
            raise OSError("Simulated mid-generation failure")

        with patch.object(G, "next_faces", interrupted):
            with self.assertRaises(OSError):
                E.run(state, 10, lambda *_: None, max_steps=1)
        self.assertEqual(state["agenda"], original)
        self.assertEqual(state["seen_inner"], original_seen)
        found = []
        E.run(state, 10, lambda f, c: found.append(c))
        self.assertTrue(state["exhausted"])
        self.assertEqual(len(found), 1)

    def test_sine_cache_preserves_integer_enclosures(self):
        rng = random.Random(119)
        for _ in range(1000):
            lo, hi = sorted(rng.sample(range(-C.S, C.S + 1), 2))
            raw = C.sqrt(1 - C.square(C.cosine(C.I(lo, hi))))
            first = C.sine(C.I(lo, hi))
            repeated = C.sine(C.I(lo, hi))
            self.assertEqual((raw.a, raw.b), (first.a, first.b))
            self.assertEqual((raw.a, raw.b), (repeated.a, repeated.b))

    def test_journal_ahead_of_checkpoint_deduplicates(self):
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stdout(
            io.StringIO()
        ):
            directory = Path(temporary)
            state = E.initial(SEED, MENU)
            E.save(state, directory / "state.pkl")
            journal = directory / "journal.jsonl"
            code = G.canonical_code(SEED + ((1, 3, 2),))
            append(journal, {"candidate_id": "old", "faces": code[1]})
            seen = {G.canonical_code(r["faces"]) for r in rows(journal)}

            def sink(faces, code):
                if code not in seen:
                    append(journal, {"faces": code[1]})
                    seen.add(code)

            state = E.load(directory / "state.pkl", MENU)
            E.run(state, 10, sink, directory / "state.pkl")
            self.assertEqual(len(rows(journal)), 1)
            self.assertTrue(E.load(directory / "state.pkl", MENU)["exhausted"])

    def test_lock_rejects_second_writer(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pass.lock"
            with lock(path):
                with self.assertRaises(RuntimeError):
                    with lock(path):
                        pass
            with lock(path):
                pass

    def test_cli_standalone_and_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            command = [
                sys.executable,
                "-S",
                str(HERE / "census.py"),
                "exhaust",
                "0",
                "1",
                "--cache-dir",
                temporary,
            ]
            first = subprocess.run(command, capture_output=True, text=True, check=True)
            second = subprocess.run(command, capture_output=True, text=True, check=True)
            self.assertIn("ENUMERATION EXHAUSTED", first.stdout)
            self.assertIn("ENUMERATION EXHAUSTED", second.stdout)
            done = json.loads((Path(temporary) / "pass_0_done.json").read_text())
            self.assertFalse(done["geometric_certification_complete"])
            self.assertEqual(done["candidates"], 0)

    def test_cli_rejects_fractional_workers_nonfinite_and_extra_numbers(self):
        for arguments in [
            ("parallel", "1.5"),
            ("all", "nan"),
            ("grow", "0", "inf"),
            ("passes", "1"),
        ]:
            with self.subTest(arguments=arguments), contextlib.redirect_stderr(
                io.StringIO()
            ):
                with self.assertRaises(SystemExit) as caught:
                    census.main(arguments)
                self.assertEqual(caught.exception.code, 2)

    def test_cache_source_mismatch_is_read_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = census.cache(temporary)
            path = directory / "configuration.json"
            path.write_text('{"version":"1.15"}\n')
            before = path.read_bytes()
            with self.assertRaises(ValueError):
                census.cache(temporary)
            self.assertEqual(path.read_bytes(), before)


class CoverAudit(unittest.TestCase):
    def test_hull_merge_preserves_all_box_points(self):
        info = G.analyze(SEED, MENU)
        edge = frozenset((0, 1))
        original = [{edge: (-C.S, -C.S // 2)}, {edge: (C.S // 2, C.S)}]
        parent = {"faces": SEED, "boxes": original}
        excluded, cover, count = C.solve(SEED, info, parent, time_cap=0, cover_limit=1)
        self.assertFalse(excluded)
        self.assertEqual(count, 0)
        self.assertEqual(cover["boxes"], [{edge: (-C.S, C.S)}])

    def test_unrelated_parent_cover_is_not_inherited(self):
        info = G.analyze(SEED, MENU)
        parent = {"faces": ((100, 101, 102),), "boxes": []}
        excluded, cover, _ = C.solve(SEED, info, parent, time_cap=0)
        self.assertFalse(excluded)
        self.assertTrue(cover["boxes"])

    def test_zero_box_budget_never_excludes(self):
        info = G.analyze(SEED, MENU)
        excluded, cover, count = C.solve(SEED, info, max_boxes=0)
        self.assertFalse(excluded)
        self.assertEqual(count, 0)
        self.assertTrue(cover["boxes"])


if __name__ == "__main__":
    unittest.main()
