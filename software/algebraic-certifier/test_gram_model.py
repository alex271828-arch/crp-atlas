"""Structural regression tests; numerical guides never establish identities."""
import unittest
from pathlib import Path
import numpy as np
import sympy as sp
from gram_model import build_gram_model, symmetry_blocks
from inputs import load_case
HERE=Path(__file__).resolve().parent

class GramModelTests(unittest.TestCase):
 def case(self,name):
  directory=HERE/'examples'/name
  return load_case(directory,spheres=directory/'spheres.json')
 def test_tetrahedron_exact_centered_gram(self):
  result=build_gram_model({'faces':[[0,1,2],[0,3,1],[0,2,3],[1,3,2]]})
  self.assertEqual(result['variables'],())
  self.assertEqual(result['gram'],sp.eye(4)/2-sp.ones(4)/8)
  self.assertEqual(result['witness'],{})
  blocks=symmetry_blocks(result)
  self.assertEqual(sum(b['dimension'] for b in blocks),3)
  self.assertTrue(all(b['numerical_rank'] is None for b in blocks))
 def test_centering_and_regular_face_constraints(self):
  result=build_gram_model(self.case('J88')); gram=result['gram']
  self.assertEqual((gram*sp.ones(gram.rows,1)).applyfunc(sp.expand),sp.zeros(gram.rows,1))
  for (i,j),d in result['known_constraints'].items():
   self.assertEqual(sp.expand(gram[i,i]+gram[j,j]-2*gram[i,j]-d),0)
  for relation in result['face_affine_relations']:
   self.assertEqual((gram*relation).applyfunc(sp.expand),sp.zeros(gram.rows,1))
  for p in result['permutations']:
   self.assertEqual((gram.extract(p,p)-gram).applyfunc(sp.expand),sp.zeros(gram.rows))
 def test_equations_independent_of_numeric_guide(self):
  case=self.case('J84'); original=build_gram_model(case)
  distorted=dict(case,positions=(np.asarray(case['positions'])*1.07+5).tolist())
  guided=build_gram_model(distorted)
  self.assertEqual(original['gram'],guided['gram'])
  self.assertEqual(original['linear_basis'],guided['linear_basis'])
  self.assertNotEqual(original['witness'],guided['witness'])
 def test_character_blocks_are_exactly_orthogonal(self):
  case=self.case('J88'); model=build_gram_model(case)
  blocks=symmetry_blocks(model,positions=case['positions'])
  self.assertEqual([b['numerical_rank'] for b in blocks],[1,1,1])
  for i,left in enumerate(blocks):
   actual=(left['basis'].T*model['gram']*left['basis']-left['gram']).applyfunc(sp.expand)
   self.assertEqual(actual,sp.zeros(left['dimension']))
   for right in blocks[i+1:]:
    cross=(left['basis'].T*model['gram']*right['basis']).applyfunc(sp.expand)
    self.assertEqual(cross,sp.zeros(left['dimension'],right['dimension']))
 def test_relabeling_preserves_counts_and_numerical_metric(self):
  case=self.case('J87'); n=len(case['positions']); p=[(3*i+2)%n for i in range(n)]
  points=[None]*n
  for i,j in enumerate(p): points[j]=case['positions'][i]
  other=dict(case,faces=[[p[v] for v in face] for face in reversed(case['faces'])],positions=points)
  first=build_gram_model(case); second=build_gram_model(other)
  b1=symmetry_blocks(first,positions=case['positions']); b2=symmetry_blocks(second,positions=points)
  self.assertEqual(len(first['variables']),len(second['variables']))
  self.assertEqual(sorted((b['dimension'],b['numerical_rank']) for b in b1),sorted((b['dimension'],b['numerical_rank']) for b in b2))
  gram=np.asarray(second['gram'].subs(second['witness']),dtype=float)
  centered=np.asarray(points)-np.mean(points,axis=0)
  self.assertLess(np.max(np.abs(gram-centered@centered.T)),1e-7)


class GramParameterizationTests(unittest.TestCase):
 def test_anchor_change_preserves_exact_family(self):
  from gram_model import anchor_parameterization
  directory=HERE/'examples'/'J88'
  case=load_case(directory,spheres=directory/'spheres.json')
  original=build_gram_model(case);symmetry_blocks(original,positions=case['positions'])
  transformed=anchor_parameterization(original)
  substituted=original['gram'].xreplace(transformed['original_parameter_expressions']).applyfunc(sp.expand)
  self.assertEqual(substituted,transformed['gram'])
  a=np.asarray(original['gram'].subs(original['witness']),dtype=float)
  b=np.asarray(transformed['gram'].subs(transformed['witness']),dtype=float)
  self.assertLess(np.max(np.abs(a-b)),1e-10)
 def test_polynomial_span_exposes_affine_constraint_without_degree_growth(self):
  from gram_model import polynomial_linear_span
  x,y=sp.symbols('x y',real=True);trace=[]
  equations,variables,assignment=polynomial_linear_span([x*x+y*y-5,x*x+y*y+x-7],[x,y],{x:2.,y:1.},trace)
  self.assertEqual(variables,[])
  self.assertEqual(assignment,{x:sp.Integer(2),y:sp.Integer(1)})
  self.assertTrue(all(sp.expand(e.subs(assignment))==0 for e in [x*x+y*y-5,x*x+y*y+x-7]))
 def test_rank_one_construction_contains_saved_metric(self):
  from gram_model import anchor_parameterization,rank_one_coordinate_system
  directory=HERE/'examples'/'J89'
  case=load_case(directory,spheres=directory/'spheres.json')
  model=build_gram_model(case);symmetry_blocks(model,positions=case['positions'])
  model=anchor_parameterization(model);system=rank_one_coordinate_system(model)
  self.assertTrue(all(not e.has(sp.Float) for e in system['equations']))
  self.assertLess(max(abs(float(e.evalf(subs=system['witness']))) for e in system['equations']),1e-7)
  for variable,expression in system['gram_parameter_expressions'].items():
   self.assertAlmostEqual(float(expression.evalf(subs=system['witness'])),model['witness'][variable],places=7)

if __name__=='__main__':unittest.main()
