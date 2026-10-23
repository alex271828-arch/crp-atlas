"""Saved partition corrections must be exact-input substitutions, not body overrides."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_body


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PartitionCorrectionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.cases = self.root/'cases'
        self.cases.mkdir()
        self.geometry = self.root/'geometry'
        self.geometry.mkdir()
        self.bi = self.geometry/'numerical.json'
        self.bi.write_text(json.dumps({'positions': [[0, 0, 0], [1, 0, 0]]}))
        self.faces = self.geometry/'input.json'
        self.faces.write_text(json.dumps({'faces': [[0, 1, 2]]}))
        self.original = self.cases/'066-Hmock.json'
        self.original_data = {
            'census_index': 66, 'display_id': 'Hmock',
            'name': 'CRP0-Hmock', 'input_file': str(self.bi),
            'input_sha256': file_hash(self.bi),
            'groups': [[0], [1, 2]], 'squared_radii': [0., 1.],
        }
        self.original.write_text(json.dumps(self.original_data))
        self.catalog = self.root/'partition_corrections'
        self.catalog.mkdir()
        self.replacement = self.catalog/self.original.name
        self.replacement_data = {
            **self.original_data,
            'groups': [[0, 1], [2]], 'squared_radii': [.5, 1.5],
        }
        self.replacement.write_text(json.dumps(self.replacement_data))
        self.manifest = {
            'format_version': 1,
            'corrections': [{
                'source_sha256': file_hash(self.original),
                'partition': self.replacement.name,
                'faces_sha256': file_hash(self.faces),
                'bi_sha256': file_hash(self.bi),
                'reason': 'Use the independently certified positive-radius partition.',
            }],
        }
        (self.catalog/'manifest.json').write_text(json.dumps(self.manifest))

    def corrected(self, path=None, data=None):
        return run_body.corrected_partition(
            path or self.original,
            data if data is not None else self.original_data,
            corrections_directory=self.catalog)

    def test_exact_saved_partition_selects_positive_replacement(self):
        path, data = self.corrected()
        self.assertEqual(path, self.replacement)
        self.assertEqual(data, self.replacement_data)
        self.assertTrue(all(radius > 0 for radius in data['squared_radii']))
        self.assertEqual(json.loads(self.original.read_text()), self.original_data)

    def test_existing_selector_aliases_resolve_to_same_corrected_body(self):
        resolver = run_body.corrected_partition

        def in_test_catalog(path, data):
            return resolver(path, data, corrections_directory=self.catalog)

        with patch('run_body.corrected_partition', side_effect=in_test_catalog):
            for selector in ('66', '066', 'Hmock', 'crp0-hmock'):
                with self.subTest(selector=selector):
                    path, data = run_body.select_body(selector, self.cases)
                    self.assertEqual(path, self.replacement)
                    self.assertEqual(data['census_index'], 66)
                    self.assertEqual(data['display_id'], 'Hmock')

    def test_custom_partition_for_same_body_is_not_overridden(self):
        custom = self.root/'custom-partition.json'
        custom_data = {**self.original_data, 'groups': [[0, 2], [1]],
                       'squared_radii': [.75, 1.75]}
        custom.write_text(json.dumps(custom_data))
        self.assertEqual(self.corrected(custom, custom_data), (custom, custom_data))

    def test_unrelated_body_is_not_overridden(self):
        other = self.root/'067-Hother.json'
        other_data = {**self.original_data, 'census_index': 67,
                      'display_id': 'Hother', 'name': 'CRP0-Hother'}
        other.write_text(json.dumps(other_data))
        self.assertEqual(self.corrected(other, other_data), (other, other_data))

    def test_changed_original_partition_is_not_overridden(self):
        changed = {**self.original_data, 'note': 'A user-edited partition'}
        self.original.write_text(json.dumps(changed))
        self.assertEqual(self.corrected(data=changed), (self.original, changed))

    def test_changed_bi_geometry_rejects_stale_correction(self):
        self.bi.write_text(json.dumps({'positions': [[0, 0, 0], [2, 0, 0]]}))
        with self.assertRaises(ValueError):
            self.corrected()

    def test_changed_face_input_rejects_stale_correction(self):
        self.faces.write_text(json.dumps({'faces': [[0, 2, 1]]}))
        with self.assertRaises(ValueError):
            self.corrected()

    def test_replacement_must_preserve_body_identity(self):
        for key, wrong in (('census_index', 67), ('display_id', 'Hother')):
            with self.subTest(key=key):
                self.replacement.write_text(json.dumps({**self.replacement_data,
                                                       key: wrong}))
                with self.assertRaises(ValueError):
                    self.corrected()

    def test_replacement_cannot_redirect_geometry_input(self):
        other_bi = self.geometry/'other-numerical.json'
        other_bi.write_bytes(self.bi.read_bytes())
        self.replacement.write_text(json.dumps({**self.replacement_data,
                                               'input_file': str(other_bi)}))
        with self.assertRaises(ValueError):
            self.corrected()

    def test_replacement_input_digest_must_agree(self):
        self.replacement.write_text(json.dumps({**self.replacement_data,
                                               'input_sha256': '0'*64}))
        with self.assertRaises(ValueError):
            self.corrected()

    def test_packaged_replacement_digest_detects_changed_partition(self):
        self.manifest['corrections'][0]['replacement_sha256'] = file_hash(self.replacement)
        (self.catalog/'manifest.json').write_text(json.dumps(self.manifest))
        self.assertEqual(self.corrected(), (self.replacement, self.replacement_data))
        self.replacement.write_text(json.dumps({**self.replacement_data,
                                               'groups': [[0, 2], [1]]}))
        with self.assertRaises(ValueError):
            self.corrected()

    def test_missing_catalog_leaves_existing_inputs_unchanged(self):
        result = run_body.corrected_partition(
            self.original, self.original_data,
            corrections_directory=self.root/'not-installed')
        self.assertEqual(result, (self.original, self.original_data))

    def test_empty_catalog_leaves_existing_inputs_unchanged(self):
        (self.catalog/'manifest.json').write_text(json.dumps({
            'format_version': 1, 'corrections': []}))
        self.assertEqual(self.corrected(), (self.original, self.original_data))


if __name__ == '__main__':
    unittest.main()
