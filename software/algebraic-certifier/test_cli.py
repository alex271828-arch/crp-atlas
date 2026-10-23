import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import sympy as sp

from check import expression


HERE = Path(__file__).resolve().parent


class CommandTests(unittest.TestCase):
    def test_certificate_parser_only_accepts_arithmetic(self):
        alpha = sp.Symbol('alpha')
        self.assertEqual(expression('(alpha**2+1)/(2*alpha)+sqrt(2)', {'alpha': alpha}),
                         (alpha**2+1)/(2*alpha)+sp.sqrt(2))
        for bad in ("__import__('os').system('true')", 'alpha.__class__', '[alpha][0]', 'open(1)', '1.5'):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                expression(bad, {'alpha': alpha})

    def test_hard_timeout_leaves_an_unresolved_result(self):
        with tempfile.TemporaryDirectory() as temp:
            run = subprocess.run([sys.executable, str(HERE/'solve.py'),
                '--case', str(HERE/'examples'/'J88'), '--spheres',
                str(HERE/'examples'/'J88'/'spheres.json'), '--out', temp,
                '--time-limit', '.05'], capture_output=True, text=True, timeout=5.)
            self.assertEqual(run.returncode, 2, run.stdout+run.stderr)
            result = json.loads((Path(temp)/'result.json').read_text())
            self.assertEqual(result['status'], 'UNRESOLVED')
            self.assertFalse(result['certified'])
            self.assertIn('wall-clock', result['reason'])
            self.assertFalse((Path(temp)/'radii.json').exists())


if __name__ == '__main__':
    unittest.main()
