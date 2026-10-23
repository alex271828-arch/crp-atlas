#!/usr/bin/env python3
"""Check a readable presentation against its saved exact radius certificate.

This checks an algebraic change of presentation, including the real embedding.
It does not independently certify the polyhedron: use check.py for that proof.
All arithmetic is exact, and displayed divisions are checked before cancellation.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys

import sympy as sp

from check import expression
from verify import ExactField


def _symbol(value, label):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', value):
        raise ValueError(f'{label} must be an arithmetic symbol')
    return sp.Symbol(value)


def _arithmetic(text, names, field=None):
    """Safely parse field arithmetic, retaining every syntactic denominator.

    The shared certificate parser is the only expression decoder. The AST walk
    first restricts this format to rational functions and checks denominators
    even when SymPy would cancel them or multiply the whole expression by zero.
    """
    if not isinstance(text, str):
        raise ValueError('An arithmetic expression must be a string')
    if field is not None:
        return _field_arithmetic(text, names, field)
    tree = ast.parse(text, mode='eval')

    def decode(node):
        return expression(ast.unparse(node), names)

    def nonzero(node):
        value = decode(node)
        if field is None:
            if not value.is_Rational:
                raise ValueError('A polynomial or endpoint has a nonconstant denominator')
            zero = value == 0
        else:
            if value.free_symbols - {field.alpha}:
                raise ValueError('A denominator contains a symbol outside the number field')
            zero = field.zero(value)
        if zero:
            raise ValueError('A displayed denominator vanishes in the selected number field')

    def visit(node):
        if isinstance(node, ast.Expression):
            visit(node.body)
        elif isinstance(node, ast.Constant) and type(node.value) is int:
            return
        elif isinstance(node, ast.Name) and node.id in names:
            return
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            visit(node.operand)
        elif isinstance(node, ast.BinOp):
            visit(node.left)
            visit(node.right)
            if isinstance(node.op, ast.Div):
                nonzero(node.right)
            elif isinstance(node.op, ast.Pow):
                exponent = decode(node.right)
                if not exponent.is_Integer:
                    raise ValueError('Presentation powers must have integer exponents')
                if exponent < 0:
                    nonzero(node.left)
            elif not isinstance(node.op, (ast.Add, ast.Sub, ast.Mult)):
                raise ValueError('Unsupported arithmetic in presentation')
        else:
            raise ValueError('Unsupported arithmetic in presentation')

    visit(tree)
    return expression(text, names)


def _field_arithmetic(text, names, field):
    """Evaluate the original AST in Q[alpha]/f, reducing each operation.

    A tuple represents a polynomial in T with already reduced field
    coefficients. This keeps generator substitution and factored Q small,
    while preserving undefined divisions that symbolic simplification hides.
    """
    modulus = field.poly
    zero = sp.Poly(0, field.alpha, domain=sp.QQ)
    one = sp.Poly(1, field.alpha, domain=sp.QQ)
    known = {}
    for name, value in names.items():
        if name == 'T' and value == sp.Symbol('T'):
            known[name] = (zero, one)
        else:
            known[name] = (sp.Poly(field.normal(value), field.alpha, domain=sp.QQ),)

    def trim(value):
        value = list(value)
        while len(value) > 1 and value[-1].is_zero:
            value.pop()
        return tuple(value)

    def add(left, right):
        return trim((left[i] if i < len(left) else zero) +
                    (right[i] if i < len(right) else zero)
                    for i in range(max(len(left), len(right))))

    def multiply(left, right):
        answer = [zero] * (len(left)+len(right)-1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                answer[i+j] = answer[i+j] + (a*b).rem(modulus)
        return trim(answer)

    def inverse(value):
        if len(value) != 1:
            raise ValueError('A denominator contains a symbol outside the number field')
        if value[0].is_zero:
            raise ValueError('A displayed denominator vanishes in the selected number field')
        return (value[0].invert(modulus),)

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return (sp.Poly(expression(str(node.value), {}), field.alpha, domain=sp.QQ),)
        if isinstance(node, ast.Name) and node.id in known:
            return known[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return tuple(-part for part in value) if isinstance(node.op, ast.USub) else value
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return add(left, right)
            if isinstance(node.op, ast.Sub):
                return add(left, tuple(-part for part in right))
            if isinstance(node.op, ast.Mult):
                return multiply(left, right)
            if isinstance(node.op, ast.Div):
                return multiply(left, inverse(right))
            if isinstance(node.op, ast.Pow):
                if len(right) != 1 or right[0].degree() > 0 or not right[0].nth(0).is_Integer:
                    raise ValueError('Presentation powers must have integer exponents')
                exponent = int(right[0].nth(0))
                if exponent < 0:
                    left, exponent = inverse(left), -exponent
                answer = (one,)
                while exponent:
                    if exponent % 2:
                        answer = multiply(answer, left)
                    exponent //= 2
                    if exponent:
                        left = multiply(left, left)
                return answer
        raise ValueError('Unsupported arithmetic in presentation')

    value = visit(ast.parse(text, mode='eval'))
    T = sp.Symbol('T')
    return sum(coefficient.as_expr()*T**i for i, coefficient in enumerate(value))


def _interval(value):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError('An isolating interval must contain two rational endpoints')
    result = [_arithmetic(item, {}) for item in value]
    if any(not endpoint.is_Rational for endpoint in result):
        raise ValueError('An isolating interval must have rational endpoints')
    return result


def _verify(saved_radii, presentation):
    if not isinstance(saved_radii, dict) or not isinstance(presentation, dict):
        raise ValueError('Saved radii and presentation must be objects')
    if type(presentation.get('format_version')) is not int or presentation['format_version'] != 1:
        raise ValueError('Unsupported presentation format version')

    alpha = _symbol(saved_radii['primitive_symbol'], 'Source primitive symbol')
    if str(alpha) == 'T':
        raise ValueError('The source primitive symbol cannot use the reserved Q variable T')
    source_names = {str(alpha): alpha}
    source_polynomial = _arithmetic(saved_radii['primitive_polynomial'], source_names)
    source_interval = _interval(saved_radii['isolating_interval'])
    field = ExactField(alpha, source_polynomial, source_interval)
    source_values = saved_radii['squared_radii']
    if not isinstance(source_values, list) or not source_values:
        raise ValueError('Source squared radii must be a nonempty list')
    source_radii = [field.normal(_arithmetic(value, source_names, field))
                    for value in source_values]

    generator = presentation['generator']
    if not isinstance(generator, dict):
        raise ValueError('The presentation generator must be an object')
    t = _symbol(generator['symbol'], 'Presentation generator')
    if str(t) == 'T' or re.fullmatch(r's[0-9]+', str(t)):
        raise ValueError('The presentation generator cannot use the reserved Q or radius symbols')
    polynomial = _arithmetic(generator['polynomial'], {str(t): t})
    interval = _interval(generator['isolating_interval'])
    target_field = ExactField(t, polynomial, interval)
    if target_field.poly.degree() != field.poly.degree():
        raise ValueError('The presentation generator does not have the source field degree')
    generator_map = field.normal(_arithmetic(generator['definition_in_source'], source_names, field))
    mapped_polynomial = _arithmetic(generator['polynomial'], {str(t): generator_map}, field)
    if not field.zero(mapped_polynomial):
        raise ValueError('The new generator does not satisfy its claimed polynomial')
    lo, hi = interval
    if lo == hi:
        if not field.zero(generator_map-lo):
            raise ValueError('The new generator does not match its rational embedding')
    elif field.sign(generator_map-lo) != 1 or field.sign(hi-generator_map) != 1:
        raise ValueError('The new generator does not match its claimed real embedding')

    definitions = presentation['definitions']
    if not isinstance(definitions, list):
        raise ValueError('Presentation definitions must be an ordered list')
    names = {str(t): generator_map}
    for definition in definitions:
        if not isinstance(definition, dict):
            raise ValueError('Every definition must be an object')
        name = str(_symbol(definition['symbol'], 'Definition symbol'))
        if name in names or name in {'T', str(alpha)} or re.fullmatch(r's[0-9]+', name):
            raise ValueError('Definition symbols must be distinct and not shadow reserved symbols')
        value = _arithmetic(definition['expression'], names, field)
        names[name] = field.normal(value)

    values = presentation['squared_radii']
    if not isinstance(values, list) or len(values) != len(source_radii):
        raise ValueError('The displayed squared radii have the wrong count')
    radius_names = {}
    for index, (value, expected) in enumerate(zip(values, source_radii)):
        actual = _arithmetic(value, names, field)
        if not field.zero(actual-expected):
            raise ValueError(f'Displayed squared radius {index} differs from the exact certificate')
        radius_names[f's{index}'] = expected

    # Factoring Q must preserve the ordered radii, including multiplicities.
    T = sp.Symbol('T')
    actual_Q = _arithmetic(presentation['Q_factored'], dict(names, **radius_names, T=T), field)
    actual_coefficients = sp.Poly(actual_Q, T).all_coeffs()
    expected_Q = _arithmetic('*'.join(f'(T-s{i})' for i in range(len(source_radii))),
                             dict(radius_names, T=T), field)
    expected_coefficients = sp.Poly(expected_Q, T).all_coeffs()
    if (len(actual_coefficients) != len(expected_coefficients) or
            any(not field.zero(actual-expected)
                for actual, expected in zip(actual_coefficients, expected_coefficients))):
        raise ValueError('Displayed Q differs from the product of the squared-radius factors')

    return {'verified': True, 'status': 'VERIFIED',
            'reason': 'Exact presentation equivalence, nonzero denominators, and real embedding verified',
            'field_degree': field.poly.degree(), 'squared_radii_count': len(source_radii)}


def verify_presentation(saved_radii: dict, presentation: dict) -> dict:
    """Return an exact equivalence report without trusting presentation metadata."""
    try:
        return _verify(saved_radii, presentation)
    except Exception as exc:
        return {'verified': False, 'status': 'UNVERIFIED',
                'reason': str(exc), 'exception': type(exc).__name__}


def worker(folder, presentation_file=None):
    try:
        folder = Path(folder).resolve()
        path = Path(presentation_file) if presentation_file else folder/'presentation.json'
        saved = json.loads((folder/'radii.json').read_text())
        presentation = json.loads(path.read_text())
        report = verify_presentation(saved, presentation)
    except Exception as exc:
        report = {'verified': False, 'status': 'UNVERIFIED',
                  'reason': str(exc), 'exception': type(exc).__name__}
    print(json.dumps(report, indent=2))
    return 0 if report.get('verified') else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory')
    parser.add_argument('--presentation', help='Presentation JSON; defaults to OUTPUT_DIRECTORY/presentation.json')
    parser.add_argument('--time-limit', type=float, default=30.)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not math.isfinite(args.time_limit) or args.time_limit <= 0:
        parser.error('Presentation verification time limit must be finite and positive')
    if args.worker:
        return worker(args.output_directory, args.presentation)
    command = [sys.executable, str(Path(__file__).resolve()), args.output_directory, '--worker']
    if args.presentation:
        command += ['--presentation', args.presentation]
    process = subprocess.Popen(command, start_new_session=True)
    try:
        return process.wait(timeout=args.time_limit)
    except (subprocess.TimeoutExpired, KeyboardInterrupt):
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
        print(json.dumps({'verified': False, 'status': 'UNVERIFIED',
                          'reason': 'Presentation verification interrupted or exceeded its time limit'}))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
