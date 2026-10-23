"""End-to-end checks of exact metric recovery, including adverse guides."""
from copy import deepcopy
from pathlib import Path
import unittest

import numpy as np
import sympy as sp

from gram_solve import solve_gram_model
from inputs import load_case
from verify import verify_solution

HERE = Path(__file__).resolve().parent


class GramRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case = load_case(HERE/'examples'/'J84',
                             spheres=HERE/'examples'/'J84'/'spheres.json')

    def test_rigid_motion_and_small_noise_do_not_supply_algebraic_coefficients(self):
        case = deepcopy(self.case)
        p = np.asarray(case['positions'], dtype=float)
        angle = .7123
        rotation = np.array([[np.cos(angle),-np.sin(angle),0],
                             [np.sin(angle),np.cos(angle),0],[0,0,1]])
        noise = np.random.default_rng(918).normal(0, 1e-10, p.shape)
        case['positions'] = (p @ rotation + [3.17, -8.41, .57] + noise).tolist()
        result = solve_gram_model(case, max_seconds=15)
        self.assertTrue(result['candidates'], result.get('reason'))
        candidate = result['candidates'][0]
        alpha = candidate['primitive_symbol']
        self.assertEqual(sp.Poly(candidate['primitive_polynomial'], alpha).degree(), 3)
        self.assertFalse(any(sp.sympify(e).atoms(sp.Float)
                             for e in candidate['squared_radii']))
        report = verify_solution(case['faces'], case['groups'], candidate['squared_radii'],
                    primitive_symbol=alpha,
                    primitive_polynomial=candidate['primitive_polynomial'],
                    isolating_interval=candidate['isolating_interval'])
        self.assertTrue(report['verified'], report)

    def test_inconsistent_one_sphere_partition_is_not_certified(self):
        case = deepcopy(self.case)
        case['groups'] = [list(range(len(case['positions'])))]
        result = solve_gram_model(case, max_seconds=15)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertFalse(result['candidates'])
        self.assertIn('partition is inconsistent', result['reason'])


if __name__ == '__main__':
    unittest.main()
