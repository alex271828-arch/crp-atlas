"""Exact, portable tests for compact algebraic presentations.

These examples use synthetic number fields rather than bodies or BI runs.
The presentation checker must establish identities independently of the
search engine and must reject undefined displayed arithmetic.
"""
import copy
import time
import unittest

import sympy as sp

from presentation_check import verify_presentation


ALPHA=sp.Symbol('alpha')
T=sp.Symbol('t')


def saved(polynomial, interval, radii):
    return {'primitive_symbol':'alpha','primitive_polynomial':str(polynomial),
            'isolating_interval':[str(x) for x in interval],
            'squared_radii':[str(x) for x in radii]}


def presentation(polynomial, interval, radii, *, definition='alpha', definitions=()):
    radii=list(map(str,radii))
    return {'format_version':1,
            'generator':{'symbol':'t','polynomial':str(polynomial),
                         'isolating_interval':[str(x) for x in interval],
                         'definition_in_source':definition,'description':'Synthetic test generator'},
            'definitions':list(definitions),'squared_radii':radii,
            'Q_factored':'*'.join(f'(T-({x}))' for x in radii),
            'score':1,'baseline_score':2,'method':'synthetic_test',
            'verified':True,'seconds':0.0}


class PresentationCheckerTests(unittest.TestCase):
    def setUp(self):
        self.source=saved(ALPHA**2-2,('7/5','3/2'),[ALPHA+2,2*ALPHA+5])
        self.identity=presentation(T**2-2,('7/5','3/2'),['t+2','2*t+5'])

    def assertVerified(self,source,result):
        report=verify_presentation(source,result)
        self.assertTrue(report['verified'],report)
        self.assertEqual(report['status'],'VERIFIED')

    def assertRejected(self,source,result):
        report=verify_presentation(source,result)
        self.assertFalse(report['verified'],report)
        self.assertEqual(report['status'],'UNVERIFIED')
        self.assertTrue(report.get('reason'))

    def test_exact_affine_generator_change_and_named_subexpressions(self):
        result=presentation(T**2-2*T-1,('12/5','5/2'),['q0+2','2*q0+5'],
                            definition='alpha+1',
                            definitions=[{'symbol':'q0','expression':'t-1'}])
        self.assertVerified(self.source,result)

    def test_negative_generator_embedding_is_valid_when_mapped_correctly(self):
        result=presentation(T**2-2,('-3/2','-7/5'),['2-t','5-2*t'],definition='-alpha')
        self.assertVerified(self.source,result)

    def test_wrong_real_embedding_is_rejected_despite_matching_polynomial(self):
        result=copy.deepcopy(self.identity)
        result['generator']['isolating_interval']=['-3/2','-7/5']
        self.assertRejected(self.source,result)

    def test_tampered_radii_and_Q_are_rejected_despite_verified_metadata(self):
        result=copy.deepcopy(self.identity)
        result['squared_radii'][0]='t+3'
        result['verified']=True
        self.assertRejected(self.source,result)
        result=copy.deepcopy(self.identity)
        result['Q_factored']='(T-(t+2))*(T-(2*t+5))+1'
        self.assertRejected(self.source,result)

    def test_named_Q_factors_refer_to_validated_radii_and_cannot_be_shadowed(self):
        result=copy.deepcopy(self.identity)
        result['Q_factored']='(T-s0)*(T-s1)'
        self.assertVerified(self.source,result)
        result['definitions']=[{'symbol':'s0','expression':'t+3'}]
        self.assertRejected(self.source,result)
        result=copy.deepcopy(self.identity)
        result['squared_radii'][0]='s0'
        self.assertRejected(self.source,result)

    def test_source_primitive_cannot_collapse_the_Q_indeterminate(self):
        source={'primitive_symbol':'T','primitive_polynomial':'T**2-2',
                'isolating_interval':['7/5','3/2'],'squared_radii':['T']}
        result=presentation(T**2-2,('7/5','3/2'),['t'],definition='T')
        result['Q_factored']='0'
        # A source-field variable named T must never turn formal Q(T)=T-s0
        # into the zero polynomial during coefficient comparison.
        self.assertRejected(source,result)

    def test_undefined_cyclic_and_shadowing_definitions_are_rejected(self):
        variants=[
            [{'symbol':'q0','expression':'q7+1'}],
            [{'symbol':'q0','expression':'q1+1'},{'symbol':'q1','expression':'q0-1'}],
            [{'symbol':'q0','expression':'t'},{'symbol':'q0','expression':'t+1'}],
            [{'symbol':'t','expression':'1'}],
            [{'symbol':'alpha','expression':'t'}],
        ]
        for definitions in variants:
            with self.subTest(definitions=definitions):
                result=copy.deepcopy(self.identity)
                result['definitions']=definitions
                self.assertRejected(self.source,result)
        result=copy.deepcopy(self.identity)
        result['squared_radii'][0]='q99'
        self.assertRejected(self.source,result)

    def test_zero_denominators_are_checked_before_symbolic_cancellation(self):
        source=saved(ALPHA**2-2,('7/5','3/2'),[1])
        for expression in ('(t**2-2)/(t**2-2)',
                           '1+0/(t**2-2)',
                           '1+0*(t**2-2)**(-1)'):
            with self.subTest(expression=expression):
                self.assertRejected(source,presentation(T**2-2,('7/5','3/2'),[expression]))
        result=presentation(T**2-2,('7/5','3/2'),['1'],
                            definitions=[{'symbol':'q0','expression':'0/(t**2-2)'}])
        self.assertRejected(source,result)
        source['squared_radii']=['(alpha**2-2)/(alpha**2-2)']
        self.assertRejected(source,presentation(T**2-2,('7/5','3/2'),['1']))

    def test_short_rational_function_equals_its_expanded_field_polynomial(self):
        modulus=ALPHA**5-ALPHA-1
        expanded=sp.invert(ALPHA+2,modulus,ALPHA)
        source=saved(modulus,('1','6/5'),[expanded+3])
        result=presentation(T**5-T-1,('1','6/5'),['3+1/(t+2)'])
        self.assertVerified(source,result)

    def test_rational_field_and_repeated_radii_preserve_Q_multiplicity(self):
        source=saved(ALPHA,('0','0'),[sp.Rational(3,8),sp.Rational(3,8)])
        result=presentation(T,('0','0'),['3/8','3/8'])
        self.assertVerified(source,result)
        result['Q_factored']='T-3/8'
        self.assertRejected(source,result)

    def test_unsafe_syntax_is_rejected_in_every_expression_position(self):
        attacks=("__import__('os').getcwd()",'t.__class__','[1][0]',
                 '(lambda: 1)()','sum([1])')
        for attack in attacks:
            for where in ('radius','definition','generator','polynomial','Q'):
                with self.subTest(attack=attack,where=where):
                    result=copy.deepcopy(self.identity)
                    if where=='radius':result['squared_radii'][0]=attack
                    elif where=='definition':result['definitions']=[{'symbol':'q0','expression':attack}]
                    elif where=='generator':result['generator']['definition_in_source']=attack
                    elif where=='polynomial':result['generator']['polynomial']=attack
                    else:result['Q_factored']=attack
                    self.assertRejected(self.source,result)


class PresentationEngineTests(unittest.TestCase):
    def engine(self):
        # Keep checker tests independently runnable while the search engine is
        # being developed; the completed suite requires this import to work.
        from presentation import build_presentation
        return build_presentation

    def test_zero_search_budget_returns_an_exact_verified_baseline(self):
        source=saved(ALPHA**3-ALPHA-1,('1','3/2'),[ALPHA**2+ALPHA+2])
        started=time.monotonic()
        result=self.engine()(source,max_seconds=0)
        self.assertLess(time.monotonic()-started,3.0)
        self.assertTrue(result['verified'],result)
        self.assertEqual(result['score'],result['baseline_score'])
        self.assertTrue(verify_presentation(source,result)['verified'])

    def test_engine_recovers_a_compact_rational_function(self):
        modulus=ALPHA**5-ALPHA-1
        expanded=sp.invert(ALPHA+2,modulus,ALPHA)
        source=saved(modulus,('1','6/5'),[expanded+3])
        checkpoints=[]
        result=self.engine()(source,max_seconds=5,
                             on_improvement=lambda value:checkpoints.append(copy.deepcopy(value)))
        self.assertTrue(verify_presentation(source,result)['verified'])
        self.assertLess(result['score'],result['baseline_score'])
        self.assertTrue(checkpoints)
        for checkpoint in checkpoints:
            self.assertTrue(verify_presentation(source,checkpoint)['verified'])
        scores=[value['score'] for value in checkpoints]
        self.assertTrue(all(a>b for a,b in zip(scores,scores[1:])))
        self.assertEqual(scores[-1],result['score'])

    def test_shared_subexpressions_are_factored_out_without_changing_radii(self):
        modulus=ALPHA**5-ALPHA-1
        shared=ALPHA**4+2*ALPHA**3+3*ALPHA**2+4*ALPHA
        source=saved(modulus,('1','6/5'),[shared+1,shared+2,shared+3,shared+4])
        result=self.engine()(source,max_seconds=5)
        self.assertTrue(verify_presentation(source,result)['verified'])
        self.assertLess(result['score'],result['baseline_score'])
        # A generator equal to the shared expression can be even shorter than
        # an auxiliary definition; either representation meets the objective.
        self.assertTrue(result['definitions'] or result['generator']['definition_in_source']!='alpha')

    def test_engine_finds_a_smaller_exact_generator_for_a_shifted_field(self):
        source=saved((ALPHA-1000)**2-2,('1001','1002'),
                     [ALPHA-998,2*ALPHA-1995])
        result=self.engine()(source,max_seconds=5)
        self.assertTrue(verify_presentation(source,result)['verified'])
        self.assertLess(result['score'],result['baseline_score'])
        self.assertNotEqual(result['generator']['definition_in_source'],'alpha')


if __name__=='__main__':unittest.main()
