"""Independent identities, geometric examples and corrupted-certificate rejection."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from algebra import Algebra, parse, ProofFailure
from reconstruct import reconstruct, cartesian
from convert import worker
from check import verify


class ConverterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temp.name)
        cls.faces = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]
        cls.case = dict(faces=cls.faces, groups=[[0, 1, 2, 3]])
        cls.saved = dict(
            primitive_symbol="alpha",
            primitive_polynomial="alpha",
            isolating_interval=["0", "0"],
            groups=cls.case["groups"],
            squared_radii=["3/8"],
            verification_method="radial",
        )
        source = cls.directory / "source"
        source.mkdir()
        out = cls.directory / "output"
        out.mkdir()
        (source / "input.json").write_text(json.dumps(cls.case))
        (source / "radii.json").write_text(json.dumps(cls.saved))
        worker(source, out)
        cls.cert = json.loads((out / "coordinates.json").read_text())

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def rejects(self, change):
        c = copy.deepcopy(self.cert)
        change(c)
        with self.assertRaises((ProofFailure, ValueError, TypeError, KeyError)):
            verify(c)

    def test_positive_constant_embedding_is_exact(self):
        A = Algebra(
            dict(
                primitive_symbol="alpha",
                primitive_polynomial="alpha**2-5",
                isolating_interval=["2", "3"],
            )
        )
        self.assertTrue(A.field.zero(A.read_base("sqrt(5)") - A.symbol))

    def test_wrong_conjugate_embedding_rejected(self):
        from unittest.mock import patch, Mock

        A = Algebra(
            dict(
                primitive_symbol="alpha",
                primitive_polynomial="alpha**2-5",
                isolating_interval=["2", "3"],
            )
        )
        with patch("algebra.sp.to_number_field", return_value=Mock(as_expr=lambda: -A.symbol)):
            with self.assertRaises(ProofFailure):
                A.read_base("sqrt(5)")

    def test_wrong_embedding_identity_rejected(self):
        from unittest.mock import patch, Mock

        A = Algebra(
            dict(
                primitive_symbol="alpha",
                primitive_polynomial="alpha**2-5",
                isolating_interval=["2", "3"],
            )
        )
        with patch("algebra.sp.to_number_field", return_value=Mock(as_expr=lambda: A.symbol + 1)):
            with self.assertRaises(ProofFailure):
                A.read_base("sqrt(5)")

    def test_nested_positive_constant_embedding(self):
        A = Algebra(
            dict(
                primitive_symbol="alpha",
                primitive_polynomial="alpha**4-10*alpha**2+5",
                isolating_interval=["3", "4"],
            )
        )
        self.assertTrue(A.field.zero(A.read_base("sqrt(5+2*sqrt(5))") - A.symbol))

    def test_tetrahedron(self):
        r = verify(self.cert)
        self.assertTrue(r["verified"])
        self.assertEqual(r["supporting_facets"], 4)

    def test_wrong_metric_coordinate(self):
        self.rejects(lambda c: c["metric_coordinates"][3].__setitem__(0, {"0": "2"}))

    def test_wrong_cartesian_coordinate(self):
        self.rejects(lambda c: c["cartesian_coordinates"][0].__setitem__(0, {}))

    def test_wrong_cartesian_map(self):
        self.rejects(lambda c: c["cartesian_transform"][0].__setitem__(0, {"0": "1"}))

    def test_wrong_center(self):
        self.rejects(lambda c: c["metric_center"].__setitem__(0, {"0": "1"}))

    def test_wrong_radius(self):
        self.rejects(lambda c: c["source_algebra"]["squared_radii"].__setitem__(0, "1"))

    def test_negative_metric(self):
        self.rejects(lambda c: c["metric"][0].__setitem__(0, "-1"))

    def test_asymmetric_metric(self):
        self.rejects(lambda c: c["metric"][0].__setitem__(1, "0"))

    def test_negative_radical(self):
        self.rejects(lambda c: c["positive_radicals"].__setitem__(0, "-1"))

    def test_unknown_radical(self):
        self.rejects(lambda c: c["metric_coordinates"][0].__setitem__(0, {"4096": "1"}))

    def test_wrong_printed_formula(self):
        self.rejects(lambda c: c["cartesian_coordinate_formulas"][0].__setitem__(0, "0"))

    def test_duplicate_sphere_label(self):
        self.rejects(lambda c: c["source_algebra"]["groups"][0].append(0))

    def test_missing_face(self):
        self.rejects(lambda c: c["input"]["faces"].pop())

    def test_float_rejected(self):
        self.rejects(lambda c: c["metric"][0].__setitem__(0, 0.375))

    def test_bool_label_rejected(self):
        self.rejects(lambda c: c["input"]["faces"][0].__setitem__(0, False))

    def test_safe_parser(self):
        for value in ['__import__("os")', "alpha.__class__", "1.0", "True", "[1]", "1/0"]:
            with self.assertRaises((ProofFailure, ValueError)):
                parse(value)

    def test_zero_denominator_before_cancellation(self):
        a = Algebra(self.saved)
        with self.assertRaises(ProofFailure):
            parse("alpha/alpha", a.symbol, a.field)

    def test_imaginary_coefficients_rejected(self):
        self.rejects(lambda c: c["metric"][0].__setitem__(1, "sqrt(-1)"))
        with self.assertRaises(ProofFailure):
            parse("sqrt(-1)**2")

    def test_radical_identities(self):
        a = Algebra(self.saved)
        x = a.sqrt(2)
        y = a.sqrt(8)
        self.assertFalse(a.sub(y, a.scale(x, 2)))
        self.assertFalse(a.sub(a.mul(x, y), a.base(4)))
        self.assertEqual(len(a.rad.generators), 1)

    def test_exact_sign(self):
        a = Algebra(self.saved)
        self.assertEqual(a.sign(a.sub(a.sqrt(2), a.base(1))), 1)
        self.assertEqual(a.sign(a.sub(a.sqrt(2), a.base(2))), -1)

    def test_reversed_orientation(self):
        case = copy.deepcopy(self.case)
        case["faces"][1].reverse()
        a, g, p, c, _ = reconstruct(case, self.saved)
        self.assertEqual(len(p), 4)

    def test_coplanar_original_faces_tile_facets(self):
        out = self.directory / "coplanar"
        out.mkdir()
        worker(ROOT / "tests/fixtures/coplanar", out)
        c = json.loads((out / "coordinates.json").read_text())
        r = verify(c)
        self.assertTrue(r["verified"])
        self.assertLess(r["supporting_facets"], r["faces"])

    def test_optimized_python_rejects_corruption(self):
        c = copy.deepcopy(self.cert)
        c["cartesian_coordinates"][0][0] = {}
        path = self.directory / "bad.json"
        path.write_text(json.dumps(c))
        p = subprocess.run(
            [sys.executable, "-O", str(ROOT / "check.py"), str(path)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(p.returncode, 2)
        self.assertIn("UNVERIFIED", p.stdout)

    def test_reducible_primitive_rejected(self):
        self.rejects(
            lambda c: c["source_algebra"].__setitem__("primitive_polynomial", "alpha**2-1")
        )

    def test_wrong_real_embedding_interval(self):
        self.rejects(lambda c: c["source_algebra"].__setitem__("isolating_interval", ["1", "2"]))

    def test_cached_sign_survives_radical_extension(self):
        a = Algebra(self.saved)
        x = a.sqrt(2)
        self.assertEqual(a.sign(x), 1)
        a.sqrt(3)
        self.assertEqual(a.sign(x), 1)

    def test_source_flags_not_used(self):
        c = copy.deepcopy(self.cert)
        c["verified"] = True
        c["status"] = "CERTIFIED"
        c["metric"][0][0] = "-1"
        with self.assertRaises(ProofFailure):
            verify(c)


if __name__ == "__main__":
    unittest.main()
