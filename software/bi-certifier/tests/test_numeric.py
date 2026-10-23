import json
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check
import metric
import run
from exact_intervals import I

FIX = Path(__file__).parent / "fixtures"


def fixture(name):
    return json.loads((FIX / name).read_text())


class CertificateTests(unittest.TestCase):
    def test_convex_metric_is_open_not_accepted(self):
        report = check.verify(fixture("tetra-root.json"))
        self.assertEqual(report["status"], "OPEN_FOR_ALGEBRAIC")
        self.assertTrue(report["convex_metric_certified"])
        self.assertFalse(report["crp_membership_certified"])
        self.assertFalse(report["exact_face_flatness_certified"])
        self.assertTrue(report["needs_algebraic_stage"])

    def test_zero_straddling_support_never_proves_convexity_or_flatness(self):
        report = check.verify(fixture("cube-root.json"))
        self.assertEqual(report["status"], "UNRESOLVED_CONVEXITY")
        self.assertTrue(report["root_certified"])
        self.assertFalse(report["excluded"])
        self.assertTrue(report["convexity"]["undecided_supports"])
        self.assertTrue(all(f["status"] == "NOT_DISPROVED" for f in report["faces"]))

    def test_nonflat_convex_root_is_excluded(self):
        report = check.verify(fixture("case-six-root.json"))
        self.assertEqual(report["status"], "EXCLUDED_NONFLAT")
        self.assertTrue(report["convex_metric_certified"])
        self.assertTrue(any(f["flatness_disproved_at_metric_root"] for f in report["faces"]))

    def test_defect_without_convexity_cannot_exclude(self):
        self.assertEqual(
            check.classify(False, [{"flatness_disproved_at_metric_root": True}]),
            "UNRESOLVED_CONVEXITY",
        )

    def test_untrusted_status_ignored(self):
        cert = fixture("cube-root.json")
        cert.update(
            status="EXCLUDED_NONFLAT",
            convex_metric_certified=True,
            exact_face_flatness_certified=True,
        )
        report = check.verify(cert)
        self.assertFalse(report["excluded"])
        self.assertFalse(report["convex_metric_certified"])

    def test_algebraic_and_symmetry_methods_not_dispatched(self):
        for method in ("exact-gram", "symmetric-affine", "orthogonal-affine"):
            with self.assertRaises(AssertionError):
                check.verify(dict(method=method))

    def test_corrupt_center_rejected(self):
        cert = fixture("tetra-root.json")
        cert["center"][0] = str(Q(cert["center"][0]) + 1)
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_zero_preconditioner_rejected(self):
        cert = fixture("tetra-root.json")
        cert["preconditioner"] = [["0"] * len(cert["center"]) for _ in cert["center"]]
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_corrupt_flip_rejected(self):
        cert = fixture("case-six-root.json")
        cert["flips"][0][2] = cert["flips"][0][0]
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_open_complex_rejected(self):
        cert = fixture("tetra-root.json")
        cert["polygons"] = cert["polygons"][:-1]
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_original_triangular_disk_chord_defect_detected(self):
        # Vertex coplanarity is automatic for three points; it does not imply
        # their distances agree with the intrinsic original triangular disk.
        report = check.measure_original_faces(
            [[0, 1, 2]], [[I(0), I(0), I(0)], [I(2), I(0), I(0)], [I(0), I(1), I(0)]]
        )
        self.assertTrue(report[0]["flatness_disproved_at_metric_root"])
        self.assertFalse(report[0]["nonzero_height_witnesses"])
        self.assertTrue(report[0]["nonzero_chord_witnesses"])

    def test_improper_input_retained(self):
        report = check.verify(fixture("improper-root.json"))
        self.assertTrue(report["root_certified"])
        self.assertEqual(report["status"], "UNRESOLVED_CONVEXITY")
        self.assertTrue(report["needs_algebraic_stage"])
        self.assertFalse(report["excluded"])

    def test_small_previous_unresolved_root_retained(self):
        report = check.verify(fixture("small-unresolved-root.json"))
        self.assertTrue(report["root_certified"])
        self.assertFalse(report["convex_metric_certified"])
        self.assertFalse(report["excluded"])


class HullMetricTests(unittest.TestCase):
    def test_cube_convex_without_flatness_claim(self):
        report = check.verify(fixture("cube-hull.json"))
        self.assertEqual(report["status"], "OPEN_FOR_ALGEBRAIC")
        self.assertTrue(report["convex_metric_certified"])
        self.assertFalse(report["exact_face_flatness_certified"])
        self.assertEqual(report["root_proof"]["patch_triangulations_checked"], [2] * 6)

    def test_improper_convex_without_merging_original_faces(self):
        cert = fixture("improper-hull.json")
        report = check.verify(cert)
        self.assertTrue(report["convex_metric_certified"])
        self.assertEqual(report["status"], "OPEN_FOR_ALGEBRAIC")
        self.assertEqual(len(report["faces"]), len(cert["polygons"]))

    def test_small_previous_unresolved_is_convex_enclosed(self):
        report = check.verify(fixture("small-unresolved-hull.json"))
        self.assertTrue(report["convex_metric_certified"])
        self.assertFalse(report["excluded"])

    def test_case_three_has_independent_interval_exclusion(self):
        report = check.verify(fixture("case-three-hull.json"))
        self.assertEqual(report["status"], "EXCLUDED_NONFLAT")
        self.assertFalse(report["internal_folds"])
        self.assertTrue(any(f["flatness_disproved_at_metric_root"] for f in report["faces"]))

    def test_missing_ridge_rejected(self):
        cert = fixture("cube-hull.json")
        cert["ridge_witnesses"].pop()
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_wrong_ridge_sign_rejected(self):
        cert = fixture("cube-hull.json")
        cert["ridge_witnesses"][0][-1] *= -1
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_lost_patch_vertex_rejected(self):
        cert = fixture("cube-hull.json")
        cert["patches"][0].pop()
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_seed_hull_wrong_connectivity_rejected(self):
        cert = fixture("cube-hull.json")
        cert["seed_triangles"][0] = cert["seed_triangles"][1]
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_wide_box_loses_neighborhood_guarantee(self):
        cert = fixture("cube-hull.json")
        cert["radius"] = "1"
        with self.assertRaises((AssertionError, ArithmeticError)):
            check.verify(cert)

    def test_hull_zero_preconditioner_rejected(self):
        cert = fixture("cube-hull.json")
        cert["preconditioner"] = [["0"] * len(cert["center"]) for _ in cert["center"]]
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_all_triangulations_enumerated_by_checker(self):
        import hull_metric as hm

        for n in range(3, 9):
            choices = hm.triangulations(tuple(range(n)))
            self.assertEqual(len(choices), hm.catalan(n - 2))
            self.assertEqual(len({frozenset(map(frozenset, ts)) for ts in choices}), len(choices))
            self.assertTrue(all(len(ts) == n - 2 for ts in choices))
        cert = fixture("cube-hull.json")
        cert["patch_triangulations_checked"] = [1] * 6  # ignored producer claim
        report = check.verify(cert)
        self.assertEqual(report["root_proof"]["patch_triangulations_checked"], [2] * 6)

    def test_piece_formulas_agree_at_flat_square(self):
        import hull_metric as hm
        from itertools import combinations

        # This tests the metric-map transition, not a coordinate template used
        # by the producer. Both diagonal choices must describe the same metric.
        p = {
            0: [hm.AD(0), hm.AD(0), hm.AD(0)],
            1: [hm.AD(1), hm.AD(0), hm.AD(0)],
            2: [hm.AD(1), hm.AD(1), hm.AD(0)],
            3: [hm.AD(0), hm.AD(1), hm.AD(0)],
        }
        edges = {metric.edge(a, b): hm.squared_distance(p[a], p[b]) for a, b in combinations(p, 2)}
        for ts in hm.triangulations((0, 1, 2, 3)):
            points = hm.develop(ts, edges)
            hm.convex_development([0, 1, 2, 3], points)
            distance = hm.squared_distance(points[0], points[2]).v
            lo, hi = map(Q, distance.exact_endpoints())
            self.assertLessEqual(lo, 2)
            self.assertGreaterEqual(hi, 2)

    def test_interval_AD_derivative_encloses_known_value(self):
        import hull_metric as hm

        x = hm.AD(I(2), {0: I(1)})
        f = (x.square() + 1).sqrt()
        derivative = f.d[0]
        # f'(2)=2/sqrt(5); compare exact rational squares.
        lo, hi = map(Q, derivative.exact_endpoints())
        self.assertLessEqual(lo * lo, Q(4, 5))
        self.assertGreaterEqual(hi * hi, Q(4, 5))


class BoundaryTests(unittest.TestCase):
    def test_metadata_coordinates_discarded(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.json"
            path.write_text(
                json.dumps(
                    dict(
                        faces=fixture("tetra-root.json")["polygons"],
                        coordinates=[[999]],
                        status="EXCLUDED_NONFLAT",
                        name="anything",
                    )
                )
            )
            item = run.load(path)[0]
            self.assertNotIn("coordinates", item)
            self.assertNotIn("status", item)

    def test_census_reads_only_candidate_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            faces = fixture("tetra-root.json")["polygons"]
            (root / "candidates_2.jsonl").write_text(
                json.dumps(dict(faces=faces, candidate_id="input-id")) + "\n"
            )
            (root / "accepted_2.jsonl").write_text("MUST NOT BE PARSED")
            (root / "ground_truth.json").write_text("MUST NOT BE PARSED")
            self.assertEqual(run.load(root, True)[0]["cid"], "input-id")

    def test_timeout_is_unresolved_and_has_handoff(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate = dict(cid="timeout-test", faces=fixture("tetra-root.json")["polygons"])
            output = Path(directory) / "case"
            with patch(
                "run.subprocess.run", side_effect=subprocess.TimeoutExpired(["metric_probe"], 1)
            ):
                result = run.attempt(candidate, output, SimpleNamespace(digits=10, timeout=1))
            self.assertEqual(result["status"], "UNRESOLVED_REALIZATION")
            self.assertFalse(result["excluded"])
            handoff = json.loads((output / "handoff.json").read_text())
            self.assertFalse(handoff["convex_metric_certified"])
            self.assertNotIn("coordinate_intervals", handoff)

    def test_proposal_timeout_salvages_only_replayed_root(self):
        actual_run = subprocess.run
        with tempfile.TemporaryDirectory() as directory:
            candidate = dict(cid="salvage", faces=fixture("cube-root.json")["polygons"])
            output = Path(directory) / "case"

            def simulate(command, **kwargs):
                if str(command[0]).endswith("metric_probe"):
                    (output / "numerical.json").write_text(
                        json.dumps([dict(positions=[[0, 0, 0]], success=False)])
                    )
                    return subprocess.CompletedProcess(command, 0)
                if len(command) > 1 and str(command[1]).endswith(("produce.py", "produce_v3.py")):
                    cert = fixture("cube-root.json")
                    cert["polygons"] = json.loads((output / "input.json").read_text())["faces"]
                    work = output / "proposal-original-metric"
                    work.mkdir()
                    (work / "root-certificate.json").write_text(json.dumps(cert))
                    raise subprocess.TimeoutExpired(command, 1)
                return actual_run(command, **kwargs)

            with patch("run.subprocess.run", side_effect=simulate):
                result = run.attempt(candidate, output, SimpleNamespace(digits=10, timeout=10))
            self.assertTrue(result["root_certified"])
            self.assertFalse(result["excluded"])
            self.assertEqual(result["status"], "UNRESOLVED_CONVEXITY")
            self.assertIn("producer_warning", result)

    def test_handoff_does_not_upgrade_root_to_convex(self):
        verified = check.verify(fixture("cube-root.json"))
        data = run.handoff(
            dict(cid="test", faces=fixture("cube-root.json")["polygons"]), verified, verified
        )
        self.assertTrue(data["root_certified"])
        self.assertFalse(data["convex_metric_certified"])
        self.assertIn("unproved", data["enclosure_target"])

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            path.write_text("preserve")
            with self.assertRaises(FileExistsError):
                run.write_json(path, {"replace": True})
            self.assertEqual(path.read_text(), "preserve")

    def test_outward_decimal(self):
        lo, hi = map(Q, run.decimals(["-1/3", "1/7"], 6)[1:-1].split(", "))
        self.assertLessEqual(lo, Q(-1, 3))
        self.assertGreaterEqual(hi, Q(1, 7))

    def test_checker_works_under_optimized_python(self):
        result = subprocess.run(
            [sys.executable, "-O", "-S", str(Path(check.__file__)), str(FIX / "tetra-root.json")],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
