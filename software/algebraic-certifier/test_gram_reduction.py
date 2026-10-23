import json
from pathlib import Path
import unittest

import sympy as sp

from gram_model import build_gram_model, symmetry_blocks
from gram_reduction import (add_reflection_patch_relations,
                            eliminate_origin_conic_pairs,
                            projected_distance_parameterization)


HERE=Path(__file__).resolve().parent


class GramReductionTests(unittest.TestCase):
    def test_line_conic_uses_known_origin_to_find_the_other_intersection(self):
        x,y,t=sp.symbols('x y t',real=True)
        equations=[x*x+y*y-2*x,y-t*x]
        trace=[]
        remaining,variables,substitution=eliminate_origin_conic_pairs(
            equations,[x,y,t],{x:1.6,y:0.8,t:0.5},trace)
        self.assertEqual(remaining,[])
        self.assertEqual(variables,[t])
        self.assertEqual(sp.cancel(substitution[x]-2/(1+t*t)),0)
        self.assertEqual(sp.cancel(substitution[y]-2*t/(1+t*t)),0)
        for equation in equations:
            self.assertEqual(sp.cancel(equation.xreplace(substitution)),0)
        self.assertTrue(trace[-1]['numerical_guide_selects_candidate_only'])

    def test_conic_without_known_origin_is_not_falsely_rationalized(self):
        x,y,t=sp.symbols('x y t',real=True)
        equations=[x*x+y*y-1,y-t*x]
        remaining,variables,substitution=eliminate_origin_conic_pairs(
            equations,[x,y,t],{x:0.8,y:0.6,t:0.75},[])
        self.assertEqual(substitution,{})
        self.assertEqual(variables,[x,y,t])
        self.assertEqual(remaining,equations)

    def test_projected_metric_anchors_come_from_exact_pair_distances(self):
        directory=HERE/'examples'/'J87'
        case={'faces':json.loads((directory/'input.json').read_text())['faces'],
              'positions':json.loads((directory/'numerical.json').read_text())['positions'],
              'groups':json.loads((directory/'spheres.json').read_text())['groups']}
        original_faces=[list(f) for f in case['faces']]
        base=build_gram_model(case)
        refined=add_reflection_patch_relations(base,case)
        self.assertLess(len(refined['variables']),len(base['variables']))
        self.assertEqual(case['faces'],original_faces)
        symmetry_blocks(refined,positions=case['positions'])
        model=projected_distance_parameterization(refined)
        anchored=[b for b in model['blocks'] if 'constant_anchor_indices' in b]
        self.assertTrue(anchored)
        for block in anchored:
            indices=block['constant_anchor_indices']
            metric=block['gram'].extract(indices,indices)
            self.assertFalse(any(v.free_symbols or v.has(sp.Float) for v in metric))
            self.assertNotEqual(metric.det(),0)
        self.assertTrue(model['protected_projection_variables'])
        self.assertFalse(any(e.has(sp.Float) for e in model['gram']))
        for (i,j),distance in base['known_constraints'].items():
            self.assertEqual(sp.expand(model['gram'][i,i]+model['gram'][j,j]
                                       -2*model['gram'][i,j]-distance),0)


if __name__=='__main__':unittest.main()
