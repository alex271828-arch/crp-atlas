"""Input integrity checks, separate from exact geometric verification."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from inputs import load_case, validate_faces


HERE = Path(__file__).resolve().parent


class InputTests(unittest.TestCase):
    def test_discovered_boundary_center_is_preserved(self):
        example = HERE/'examples'/'J92'
        sphere = json.loads((example/'spheres.json').read_text())
        with patch('numeric_partition.solve_min_spheres_numeric', return_value=sphere):
            case = load_case(example)
        self.assertEqual(case['center'], sphere['center'])
        self.assertEqual(len(case['groups']), 2)
        self.assertLess(abs(case['preflight']['center_clearance']), 1e-8)
        self.assertTrue(any('retaining its partition' in note for note in case['notes']))

    def test_worked_examples_have_interior_centers(self):
        for name in ('J84', 'J88', 'CRP87', 'CRP140'):
            with self.subTest(name=name):
                path = HERE/'examples'/name
                case = load_case(path, spheres=path/'spheres.json')
                self.assertGreater(case['preflight']['center_clearance'], .4)
                self.assertEqual(case['preflight']['within_shell_center_rank'], 3)

    def test_closed_sphere_required(self):
        tetra = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]]
        self.assertEqual(validate_faces(tetra), tetra)
        with self.assertRaisesRegex(ValueError, 'two faces'):
            validate_faces(tetra[:-1])
        disjoint = tetra + [[v+4 for v in f] for f in tetra]
        with self.assertRaisesRegex(ValueError, 'Euler'):
            validate_faces(disjoint)

    def test_wrong_sphere_membership_is_not_silently_used(self):
        example = HERE/'examples'/'J84'
        sphere = json.loads((example/'spheres.json').read_text())
        a, b = sphere['groups']
        a[0], b[0] = b[0], a[0]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'spheres.json'
            path.write_text(json.dumps(sphere))
            with self.assertRaisesRegex(ValueError, 'does not match'):
                load_case(example, spheres=path)

    def test_sphere_source_hash_is_checked(self):
        example = HERE/'examples'/'J84'
        sphere = json.loads((example/'spheres.json').read_text())
        sphere['input_sha256'] = '0'*64
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'spheres.json'
            path.write_text(json.dumps(sphere))
            with self.assertRaisesRegex(ValueError, 'SHA256'):
                load_case(example, spheres=path)


if __name__ == '__main__':
    unittest.main()
