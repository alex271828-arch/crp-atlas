"""Behavioral checks for symmetry reduction, gradient descent, and failure status."""

import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import metric
import run
import symmetry


def johnson(n):
    return json.loads((ROOT / "examples/johnson92.json").read_text())[n - 1]["faces"]


def example(name):
    return json.loads((ROOT / "examples" / name).read_text())["faces"]


def numerical(faces, limit=50000, starts=3, program="metric_probe"):
    faces = metric.orient_polygons(faces)
    group = symmetry.compute(faces)
    with tempfile.TemporaryDirectory() as tmp:
        source, output = Path(tmp) / "metric.txt", Path(tmp) / "numerical.json"
        symmetry.write_metric(source, faces, group["vertex_orbits"])
        subprocess.run(
            [
                str(ROOT / program),
                str(source),
                str(output),
                str(limit),
                "1e-9",
                "1.5",
                *([str(starts)] if program == "metric_probe" else []),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=20,
        )
        return json.loads(output.read_text())[0]


class SymmetryTests(unittest.TestCase):
    def test_full_groups_include_reflections(self):
        for name, order in [("tetrahedron.json", 24), ("cube.json", 48), ("octahedron.json", 48)]:
            with self.subTest(name=name):
                group = symmetry.compute(example(name))
                self.assertEqual(group["order"], order)
                self.assertEqual(len(group["vertex_orbits"]), 1)

    def test_square_pyramid_orbits(self):
        group = symmetry.compute(johnson(1))
        self.assertEqual(group["order"], 8)
        self.assertEqual(sorted(map(len, group["vertex_orbits"])), [1, 4])

    def test_relabeling_and_cycle_start_do_not_change_orbits(self):
        faces = johnson(88)
        original = symmetry.compute(faces)
        labels = list(range(sum(map(len, original["vertex_orbits"]))))
        random.Random(881).shuffle(labels)
        transformed = [[labels[v] for v in f[1:] + f[:1]] for f in reversed(faces)]
        changed = symmetry.compute(transformed)
        self.assertEqual(original["order"], changed["order"])
        expected = sorted(sorted(labels[v] for v in orbit) for orbit in original["vertex_orbits"])
        self.assertEqual(expected, changed["vertex_orbits"])


@unittest.skipUnless((ROOT / "metric_probe").is_file(), "Build the numerical helper first")
class GradientTests(unittest.TestCase):
    def test_known_tetrahedron_radius_and_curvature(self):
        data = numerical(example("tetrahedron.json"))
        self.assertTrue(data["success"])
        self.assertLessEqual(data["max_radial_defect"], 1e-9)
        self.assertAlmostEqual(data["radii"][0], (3 / 8) ** 0.5, places=9)

    def test_orbit_equality_and_monotone_descent(self):
        for n in (5, 6, 27, 30, 48, 64, 76, 84, 88):
            with self.subTest(johnson=n):
                data = numerical(johnson(n))
                self.assertTrue(data["success"], data["reason"])
                for orbit in data["vertex_orbits"]:
                    self.assertEqual(len({data["radii"][v] for v in orbit}), 1)
                    self.assertEqual(len({data["squared_radii"][v] for v in orbit}), 1)
                self.assertEqual(data["max_accepted_energy_increase"], 0)
                energies = [sample[1] for sample in data["energy_trace"]]
                self.assertTrue(all(b < a for a, b in zip(energies, energies[1:])))
                self.assertLess(data["curvature_energy"], data["initial_curvature_energy"])
                self.assertLess(data["max_original_face_pair_error"], 1e-7)

    def test_stall_stops_and_equal_radius_restart_recovers(self):
        single = numerical(johnson(30), starts=1)
        self.assertFalse(single["success"])
        self.assertIn("stalled", single["reason"])
        self.assertLess(single["iterations"], 100)
        recovered = numerical(johnson(30))
        self.assertTrue(recovered["success"])
        self.assertEqual(recovered["selected_attempt"], 1)
        self.assertEqual(len(recovered["attempts"]), 2)

    def test_exhausted_budget_never_claims_convergence(self):
        data = numerical(example("tetrahedron.json"), limit=0, starts=1)
        self.assertFalse(data["success"])
        self.assertEqual(data["iterations"], 0)
        self.assertIn("budget exhausted", data["reason"])

    def test_numerical_only_does_not_certify(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "case"
            args = SimpleNamespace(
                digits=10, timeout=30, numerical_only=True, max_iterations=0, radius_starts=1
            )
            result = run.attempt(
                dict(cid="tetra", faces=example("tetrahedron.json")), directory, args
            )
            self.assertEqual(result["status"], "UNRESOLVED_NUMERICAL")
            self.assertFalse(result["excluded"])
            self.assertFalse(result["root_certified"])
            self.assertFalse(result["convex_metric_certified"])
            self.assertFalse((directory / "proposal.log").exists())
            self.assertTrue((directory / "symmetry.json").exists())

    def test_original_metric_proposal_avoids_parallel_trace(self):
        import contextlib
        import io
        import check
        import produce

        faces = metric.orient_polygons(johnson(7))
        data = numerical(faces)
        with self.assertRaisesRegex(AssertionError, "parallel intrinsic edges"):
            metric.replay_flips(faces, data["flips"])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "input.json").write_text(json.dumps(dict(faces=faces)))
            original = json.dumps([data])
            (directory / "numerical.json").write_text(original)
            with contextlib.redirect_stdout(io.StringIO()):
                produce.produce(directory, 10)
            chosen = (directory / "chosen-certificate.txt").read_text().strip()
            certificate = json.loads((directory / chosen).read_text())
            self.assertEqual(certificate["flips"], [])
            self.assertTrue(check.verify(certificate)["convex_metric_certified"])
            self.assertEqual((directory / "numerical.json").read_text(), original)

    @unittest.skipUnless((ROOT / "gradient_check").is_file(), "Build the gradient diagnostic first")
    def test_analytic_gradient_against_independent_difference(self):
        for n in (1, 5, 6, 27, 30, 76, 84, 88):
            with self.subTest(johnson=n):
                data = numerical(johnson(n), program="gradient_check")
                self.assertLess(data["gradient_relative_error"], 1e-7)
                self.assertLess(data["replay_error"], 1e-12)


if __name__ == "__main__":
    unittest.main()
