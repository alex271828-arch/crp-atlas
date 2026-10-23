"""Saved arbitrary-center certificates use safe parsing and exact rechecking."""
import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from check import geometric_expressions, worker


FACES = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]]
GEOMETRY = {'metric': [['1', '1/2', '1/2'], ['1/2', '1', '1/2'], ['1/2', '1/2', '1']],
            'coordinates': [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'center': [1, 1, 1]}


class SavedGeometryTests(unittest.TestCase):
    @staticmethod
    def saved():
        return {'primitive_symbol': 'alpha', 'primitive_polynomial': 'alpha',
                'isolating_interval': ['0', '0'], 'groups': [[0], [1, 2, 3]],
                'squared_radii': ['6', '3'], 'Q_coefficients': ['1', '-9', '18'],
                'verification_method': 'exact_gram_geometry',
                'geometric_certificate': copy.deepcopy(GEOMETRY)}

    def recheck(self, saved):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'input.json').write_text(json.dumps({'faces': FACES, 'groups': saved['groups']}))
            (folder/'radii.json').write_text(json.dumps(saved))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                status = worker(folder)
            return status, json.loads(output.getvalue())

    def test_exterior_center_certificate_is_rechecked_exactly(self):
        status, report = self.recheck(self.saved())
        self.assertEqual(status, 0, report)
        self.assertTrue(report['certified'])
        self.assertEqual(report['verification_method'], 'exact_gram_geometry')
        self.assertEqual(report['center_location'], 'exterior')
        self.assertEqual(report['checks']['exact_sphere_memberships'], 4)
        self.assertEqual(report['checks']['strictly_supporting_face_planes'], 4)

    def test_old_radial_certificate_remains_supported(self):
        saved = self.saved()
        del saved['verification_method']
        del saved['geometric_certificate']
        saved.update(groups=[[0, 1, 2, 3]], squared_radii=['3/8'],
                     Q_coefficients=['1', '-3/8'])
        status, report = self.recheck(saved)
        self.assertEqual(status, 0, report)
        self.assertTrue(report['certified'])
        self.assertEqual(report['verification_method'], 'radial')

    def test_changed_coordinates_or_center_do_not_certify(self):
        for key, value in [('coordinates', [[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 1]]),
                           ('center', [2, 1, 1])]:
            with self.subTest(key=key):
                saved = self.saved()
                saved['geometric_certificate'][key] = value
                status, report = self.recheck(saved)
                self.assertEqual(status, 2)
                self.assertFalse(report['certified'])

    def test_q_coefficients_are_checked_for_coordinate_certificates(self):
        saved = self.saved()
        saved['Q_coefficients'][-1] = '17'
        status, report = self.recheck(saved)
        self.assertEqual(status, 2)
        self.assertFalse(report['certified'])
        self.assertIn('Q coefficients', report['reason'])

    def test_geometry_expression_strings_cannot_execute_code(self):
        saved = self.saved()
        saved['geometric_certificate']['metric'][0][0] = "__import__('os').system('echo unsafe')"
        with patch('os.system', side_effect=AssertionError('must not execute')) as command:
            status, report = self.recheck(saved)
        self.assertEqual(status, 2)
        self.assertFalse(report['certified'])
        self.assertIn('Unsupported expression syntax', report['reason'])
        command.assert_not_called()

    def test_malformed_geometry_and_unknown_method_are_rejected(self):
        for replacement in ({}, {'metric': [], 'coordinates': [], 'center': []}, None):
            saved = self.saved()
            saved['geometric_certificate'] = replacement
            status, report = self.recheck(saved)
            self.assertEqual(status, 2)
            self.assertFalse(report['certified'])
        saved = self.saved()
        saved['verification_method'] = 'trust_saved_claim'
        status, report = self.recheck(saved)
        self.assertEqual(status, 2)
        self.assertFalse(report['certified'])
        self.assertIn('Unsupported saved verification method', report['reason'])

    def test_float_entries_are_rejected_by_the_arithmetic_reader(self):
        payload = copy.deepcopy(GEOMETRY)
        payload['center'][0] = 1.0
        with self.assertRaisesRegex(ValueError, 'Unsupported expression syntax'):
            geometric_expressions(payload, {})


if __name__ == '__main__':
    unittest.main()
