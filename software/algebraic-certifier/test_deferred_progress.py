"""Deferred angle closures retain exact progress between reduction passes."""
import unittest
from unittest.mock import patch

import sympy as sp

import encoder
import reduce


def unit_corner(cosine):
    return dict(N=cosine, M=sp.S.One, S=1-cosine**2,
                D_left=sp.S.One, D_right=sp.S.One)


class DeferredProgressTests(unittest.TestCase):
    def setUp(self):
        self.budget = reduce._time_limit(15)
        self.budget.__enter__()
        self.addCleanup(self.budget.__exit__, None, None, None)
        self.x, self.y, self.z = sp.symbols('x y z', real=True)
        self.radii = (self.x, self.y, self.z)
        # The five-angle fan is initially too large to materialize. The
        # equality of the first and third angles gives z=x+2 cheaply.
        self.fan = dict(counts=(3, 1, 1), rhs_pi=sp.Integer(2))
        self.equal_angles = dict(counts=(1, 0, -1), rhs_pi=sp.S.Zero)

    def model(self, cosine, relations, short=()):
        return dict(radius_variables=self.radii, variables=self.radii,
                    witness={v: float(offset+cosine)
                             for v, offset in zip(self.radii, (2, 3, 4))},
                    consequence_equations=[self.y-self.x-1],
                    corner_types=[dict(unit_corner(v-offset),
                                       y=sp.Symbol('corner_sine_'+str(offset), positive=True))
                                  for v, offset in zip(self.radii, (2, 3, 4))],
                    deferred_relations=relations,
                    short_radical_equations=list(short))

    def assert_candidate(self, result, cosine):
        self.assertEqual(result['status'], 'ALGEBRAIC_CANDIDATE', result.get('reason'))
        values = result['candidate']['algebraic_values']
        for variable, offset in zip(self.radii, (2, 3, 4)):
            self.assertEqual(sp.to_number_field(values[variable]-offset-cosine).coeffs(), [0])

    def test_later_fan_reduction_unlocks_earlier_fan_on_next_pass(self):
        cosine = (sp.sqrt(5)-1)/4
        model = self.model(cosine, [self.fan, self.equal_angles])
        # Keep this a genuine multi-pass case under the new bivariate
        # allowance: the edge product is not visibly paired until z=x+2
        # turns this positive scale into 1. The represented cosine is x-2.
        scale = self.z-self.x-1
        value = self.x-2
        model['corner_types'][0].update(
            N=scale*value, M=scale**2, S=scale**2*(1-value**2),
            D_left=scale**2, D_right=sp.S.One)
        self.assertEqual(encoder.materialize_relation(
            self.fan, model['corner_types'])[0], [])
        # Use real angle identities and exact elimination throughout. Count
        # translations only to detect regenerating an already inherited fan.
        with patch('encoder.materialize_relation', wraps=encoder.materialize_relation) as materialize:
            result = reduce.solve_model(model, max_seconds=10)
        self.assert_candidate(result, cosine)
        self.assertEqual(sum(call.args[0] == self.equal_angles
                             for call in materialize.call_args_list), 1)
        instantiations = [row for row in result['derivation']
                          if row['operation'] == 'local_angle_template_instantiation'
                          and row['counts'] == self.fan['counts']]
        self.assertGreater(len(instantiations), 1)
        self.assertTrue(instantiations[-1]['substituted_radius_family_first'])
        self.assertFalse(any(phase['name'] == 'compact_signed_radical_relations'
                             for phase in result['phases']))

    def test_partial_deferred_reduction_is_available_to_compact_fallback(self):
        cosine = sp.Rational(1, 2)
        # The deferred phase cannot isolate x. The compact equation can do
        # so only if it inherits the other two radii's derived relations.
        model = self.model(cosine, [self.equal_angles], [self.x-sp.Rational(5, 2)])
        result = reduce.solve_model(model, max_seconds=10)
        self.assert_candidate(result, cosine)
        closure = next(phase for phase in result['phases']
                       if phase['name'] == 'deferred_local_angle_closure')
        self.assertEqual(closure['solutions'], 0)
        self.assertTrue(any(phase['name'] == 'compact_signed_radical_relations'
                            and phase['solutions'] > 0 for phase in result['phases']))

    def test_unsupported_fan_stops_when_no_parameter_is_eliminated(self):
        model = self.model((sp.sqrt(5)-1)/4, [self.fan])
        with patch('encoder.materialize_relation', wraps=encoder.materialize_relation) as materialize:
            result = reduce.solve_model(model, max_seconds=10)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertFalse(result.get('timed_out', False), result.get('reason'))
        # One original template and one bounded two-parameter retry; no
        # subsequent pass is attempted without a decrease in dimension.
        self.assertEqual(materialize.call_count, 2)
        self.assertFalse(any(row['operation'] == 'phase_budget_exhausted'
                             for row in result['derivation']))

    def test_timeout_after_completed_reduction_keeps_it_for_fallback(self):
        model = self.model(sp.Rational(1, 2), [self.equal_angles],
                           [self.x-sp.Rational(5, 2)])
        original = reduce._groebner_solutions
        calls = 0

        def bounded_closure(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise reduce.ReductionTimeout('synthetic timeout in deferred closure')
            return original(*args, **kwargs)

        with patch('reduce._groebner_solutions', side_effect=bounded_closure):
            result = reduce.solve_model(model, max_seconds=10)
        self.assert_candidate(result, sp.Rational(1, 2))
        self.assertTrue(any(row['operation'] == 'phase_budget_exhausted'
                            and row['phase'] == 'deferred_local_angle_closure'
                            for row in result['derivation']))

    def test_timeout_during_reduction_does_not_commit_incomplete_pass(self):
        model = self.model(sp.Rational(1, 2), [self.equal_angles],
                           [self.x-sp.Rational(5, 2)])
        original = reduce._staged_linear
        calls = 0

        def interrupted_reduction(*args, **kwargs):
            nonlocal calls
            calls += 1
            reduced = original(*args, **kwargs)
            if calls == 2:
                # Even if individual pivots reached the trace, the pass
                # did not return a complete and consistent state.
                raise reduce.ReductionTimeout('synthetic interruption before commit')
            return reduced

        with patch('reduce._staged_linear', side_effect=interrupted_reduction), \
                patch('reduce._radical_square_classes', wraps=reduce._radical_square_classes) as radicals:
            result = reduce.solve_model(model, max_seconds=10)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertTrue(any(row['operation'] == 'phase_budget_exhausted'
                            and row['phase'] == 'deferred_local_angle_closure'
                            for row in result['derivation']))
        # Fallback must see the last completed state: y=x+1, while the
        # uncommitted z=x+2 relation has left z dependent on a parameter.
        assignments = radicals.call_args.args[1]
        self.assertEqual(sp.cancel(assignments[self.y]-assignments[self.x]-1), 0)
        self.assertTrue(sp.cancel(assignments[self.z]-assignments[self.x]-2).free_symbols)


if __name__ == '__main__':
    unittest.main()
