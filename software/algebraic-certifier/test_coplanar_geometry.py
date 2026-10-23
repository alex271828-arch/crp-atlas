"""Coplanar regular faces must certify as exact subdivisions of hull facets."""
import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile
import unittest

import sympy as sp

from check import worker
from verify_gram import verify_gram_solution


BODY37_FACES = [[5, 4, 6], [0, 3, 4, 1], [7, 2, 6], [0, 2, 3],
                [4, 3, 7], [5, 1, 4], [0, 1, 2], [3, 2, 7],
                [2, 1, 5, 6], [6, 4, 7]]
BODY37_GROUPS = [[0, 1, 2, 3, 4], [7], [5, 6]]
BODY37_RADII = [sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(11, 6)]
BODY37_GEOMETRY = {
    'metric': sp.diag(1, 1, sp.Rational(1, 2)).tolist(),
    'coordinates': [[0, 0, 0], [0, 1, 0], [sp.Rational(1, 2), sp.Rational(1, 2), 1],
                    [1, 0, 0], [1, 1, 0],
                    [sp.Rational(1, 2), sp.Rational(11, 6), sp.Rational(1, 3)],
                    [1, sp.Rational(4, 3), sp.Rational(4, 3)],
                    [sp.Rational(3, 2), sp.Rational(1, 2), 1]],
    'center': [sp.Rational(1, 2), sp.Rational(1, 2), 0],
}


class CoplanarGeometryTests(unittest.TestCase):
    def assertBody37Certified(self, report):
        self.assertTrue(report['verified'], report)
        self.assertTrue(report['certified'], report)
        self.assertEqual(report['status'], 'CERTIFIED')
        self.assertEqual(report['center_location'], 'boundary')
        self.assertEqual(report['checks']['exact_sphere_memberships'], 8)
        self.assertEqual(report['checks']['coplanar_subdivided_facets'], 2)
        self.assertEqual(report['checks']['maximal_supporting_facets'], 8)
        self.assertEqual(report['checks']['supporting_face_planes'], 10)
        self.assertEqual(report['checks']['strictly_supporting_face_planes'], 6)
        self.assertTrue(report['checks']['embedded_closed_convex_boundary'])

    def assertRejected(self, report):
        self.assertFalse(report['verified'], report)
        self.assertFalse(report['certified'], report)
        self.assertEqual(report['status'], 'UNVERIFIED')
        self.assertIn('reason', report)

    def verify(self, faces=None, radii=None, geometry=None):
        return verify_gram_solution(
            BODY37_FACES if faces is None else faces,
            BODY37_GROUPS,
            BODY37_RADII if radii is None else radii,
            BODY37_GEOMETRY if geometry is None else geometry,
        )

    def test_body37_exact_coplanar_subdivision_is_certified(self):
        self.assertBody37Certified(self.verify())

    def test_independently_reversing_each_face_preserves_certification(self):
        for face_id in range(len(BODY37_FACES)):
            with self.subTest(face=face_id):
                faces = copy.deepcopy(BODY37_FACES)
                faces[face_id].reverse()
                self.assertBody37Certified(self.verify(faces=faces))

    def test_face_order_and_mixed_orientation_preserve_certification(self):
        faces = [list(reversed(face)) if index % 2 else list(face)
                 for index, face in enumerate(reversed(BODY37_FACES))]
        self.assertBody37Certified(self.verify(faces=faces))

    def test_changed_squared_radii_do_not_certify(self):
        for group_id in range(len(BODY37_RADII)):
            with self.subTest(group=group_id):
                radii = list(BODY37_RADII)
                radii[group_id] += sp.Rational(1, 6)
                self.assertRejected(self.verify(radii=radii))

    def test_changed_coordinates_do_not_certify(self):
        for vertex in (2, 5):
            with self.subTest(vertex=vertex):
                geometry = copy.deepcopy(BODY37_GEOMETRY)
                geometry['coordinates'][vertex][2] += sp.Rational(1, 6)
                self.assertRejected(self.verify(geometry=geometry))

    @staticmethod
    def saved_certificate():
        return {
            'primitive_symbol': 'alpha',
            'primitive_polynomial': 'alpha',
            'isolating_interval': ['0', '0'],
            'groups': copy.deepcopy(BODY37_GROUPS),
            'squared_radii': [str(radius) for radius in BODY37_RADII],
            'Q_coefficients': ['1', '-23/6', '53/12', '-11/8'],
            'verification_method': 'exact_gram_geometry',
            'geometric_certificate': copy.deepcopy(BODY37_GEOMETRY),
        }

    def recheck(self, saved):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder / 'input.json').write_text(json.dumps({
                'faces': BODY37_FACES, 'groups': BODY37_GROUPS,
            }))
            (folder / 'radii.json').write_text(json.dumps(saved, default=str))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                status = worker(folder)
            return status, json.loads(output.getvalue())

    def test_saved_coplanar_certificate_round_trip(self):
        status, report = self.recheck(self.saved_certificate())
        self.assertEqual(status, 0, report)
        self.assertBody37Certified(report)
        self.assertEqual(report['verification_method'], 'exact_gram_geometry')

    def test_saved_coplanar_certificate_recomputes_geometry_and_radii(self):
        for mutation in ('coordinate', 'radius', 'polynomial'):
            with self.subTest(mutation=mutation):
                saved = self.saved_certificate()
                if mutation == 'coordinate':
                    saved['geometric_certificate']['coordinates'][2][2] = '7/6'
                elif mutation == 'radius':
                    saved['squared_radii'][0] = '2/3'
                else:
                    saved['Q_coefficients'][-1] = '-1'
                status, report = self.recheck(saved)
                self.assertEqual(status, 2, report)
                self.assertRejected(report)


if __name__ == '__main__':
    unittest.main()
