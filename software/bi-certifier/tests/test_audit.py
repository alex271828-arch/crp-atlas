"""Independent arithmetic and adversarial checks for the version 4 audit."""

from fractions import Fraction as Q
import json
from pathlib import Path
import random
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check
import exact_intervals as iv
from field import chord_squared
from proof import precision

FIX = Path(__file__).parent / "fixtures"


def fixture(name="tetra-root.json"):
    return json.loads((FIX / name).read_text())


class ArithmeticTests(unittest.TestCase):
    def tearDown(self):
        precision(160)

    def test_outward_operations_against_rational_endpoints(self):
        rng = random.Random(4426)
        for bits in (80, 160, 257):
            precision(bits)
            for _ in range(100):
                a, b = sorted(Q(rng.randrange(-50, 51), rng.randrange(1, 20)) for _ in range(2))
                c, d = sorted(Q(rng.randrange(-50, 51), rng.randrange(1, 20)) for _ in range(2))
                x, y = iv.I(a, b), iv.I(c, d)
                products = [u * v for u in (a, b) for v in (c, d)]
                for result, lo, hi in [
                    (x + y, a + c, b + d),
                    (x - y, a - d, b - c),
                    (x * y, min(products), max(products)),
                ]:
                    actual_lo, actual_hi = map(Q, result.exact_endpoints())
                    self.assertLessEqual(actual_lo, lo)
                    self.assertGreaterEqual(actual_hi, hi)
                if c > 0 or d < 0:
                    quotients = [u / v for u in (a, b) for v in (c, d)]
                    lo, hi = map(Q, (x / y).exact_endpoints())
                    self.assertLessEqual(lo, min(quotients))
                    self.assertGreaterEqual(hi, max(quotients))
                z = x.square()
                lo, hi = map(Q, z.sqrt().exact_endpoints())
                lower = 0 if a <= 0 <= b else min(abs(a), abs(b))
                self.assertLessEqual(lo, lower)
                self.assertGreaterEqual(hi, max(abs(a), abs(b)))

    def test_integer_endpoints_are_exact_at_every_precision(self):
        for bits in (80, 160, 4096):
            precision(bits)
            for a, b in ((0, 0), (-7, 13), (2**250, 2**250 + 1)):
                interval = iv.I(a, b)
                self.assertEqual(interval.lo, a * iv.SCALE)
                self.assertEqual(interval.hi, b * iv.SCALE)
            with self.assertRaises(ValueError):
                iv.I(1, -1)

    def test_dot_product_matches_exact_endpoint_extrema(self):
        rng = random.Random(182)
        for n in (0, 1, 3, 37):
            pairs = []
            for _ in range(n):
                a, b = sorted(rng.randrange(-1000, 1000) for _ in range(2))
                c, d = sorted(rng.randrange(-1000, 1000) for _ in range(2))
                pairs.append((iv.I(a, b, raw=True), iv.I(c, d, raw=True)))
            low = high = Q(0)
            for a, b in pairs:
                products = [Q(x * y, iv.SCALE**2) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
                low += min(products)
                high += max(products)
            self.assertEqual(
                iv.sum_products(pairs).exact_endpoints(), iv.I(low, high).exact_endpoints()
            )
            old = sum((a * b for a, b in pairs), iv.I(0))
            new = iv.sum_products(pairs)
            self.assertGreaterEqual(new.lo, old.lo)
            self.assertLessEqual(new.hi, old.hi)

    def test_exact_constant_cache_survives_precision_changes(self):
        precision(80)
        a = chord_squared(5, 2).interval()
        alo, ahi = map(Q, a.exact_endpoints())
        precision(257)
        b = chord_squared(5, 2).interval()
        blo, bhi = map(Q, b.exact_endpoints())
        self.assertGreaterEqual(blo, alo)
        self.assertLessEqual(bhi, ahi)
        self.assertLess(bhi - blo, ahi - alo)

    def test_reversed_sub_dyadic_interval_rejected(self):
        with self.assertRaises(ValueError):
            iv.I(Q(1, 2**200), Q(0))
        with self.assertRaises(TypeError):
            iv.I(0.1, raw=True)


class AdversarialTests(unittest.TestCase):
    def test_rational_certificate_inputs_reject_float_and_boolean(self):
        for field in ("center", "radius", "preconditioner"):
            for value in (0.0, True):
                cert = fixture()
                if field == "center":
                    cert[field][0] = value
                elif field == "radius":
                    cert[field] = value
                else:
                    cert[field][0][0] = value
                with self.subTest(field=field, value=value), self.assertRaises(AssertionError):
                    check.verify(cert)

    def test_duplicate_and_bad_triangle_labels_rejected(self):
        for key in ("triangles", "seed_triangles"):
            for mutation in ("duplicate", "float", "negative", "outside"):
                cert = fixture("cube-hull.json")
                if mutation == "duplicate":
                    cert[key].append(cert[key][0])
                else:
                    cert[key][0][0] = {"float": 0.0, "negative": -1, "outside": 1000}[mutation]
                with self.subTest(key=key, mutation=mutation), self.assertRaises(AssertionError):
                    check.verify(cert)

    def test_invalid_base_and_precision_rejected(self):
        for base in ([0, 0, 1], [False, 1, 2], [0.0, 1, 2], [-1, 1, 2], [0, 1]):
            cert = fixture()
            cert["base"] = base
            with self.subTest(base=base), self.assertRaises(AssertionError):
                check.verify(cert)
        for bits in (True, 79, 4097, 160.0):
            cert = fixture()
            cert["bits"] = bits
            with self.subTest(bits=bits), self.assertRaises(AssertionError):
                check.verify(cert)

    def test_input_proof_flags_never_override_hull_checks(self):
        cert = fixture("cube-hull.json")
        cert.update(verified=True, status="OPEN_FOR_ALGEBRAIC", convex_metric_certified=True)
        cert["ridge_witnesses"][0][-1] *= -1
        with self.assertRaises(AssertionError):
            check.verify(cert)

    def test_optimized_python_still_rejects_invalid_proof(self):
        import tempfile

        with tempfile.TemporaryDirectory() as temp:
            cert = fixture()
            cert["center"][0] = "100"
            path = Path(temp) / "bad.json"
            path.write_text(json.dumps(cert))
            result = subprocess.run(
                [sys.executable, "-O", "-S", str(Path(check.__file__)), str(path)],
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b"CertificateError", result.stderr)


if __name__ == "__main__":
    unittest.main()
