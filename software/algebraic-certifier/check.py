#!/usr/bin/env python3
"""Independently recheck saved exact certificates; no numerical data are used."""
import argparse
import ast
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import sympy as sp


def expression(text, names):
    """Read arithmetic expressions without eval, sympify(str), or parse_expr."""
    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return names[node.id]
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return -visit(node.operand)
            if isinstance(node.op, ast.UAdd):
                return visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return a+b
            if isinstance(node.op, ast.Sub):
                return a-b
            if isinstance(node.op, ast.Mult):
                return a*b
            if isinstance(node.op, ast.Div):
                return a/b
            if isinstance(node.op, ast.Pow) and b.is_Rational:
                return a**b
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == 'sqrt' and len(node.args) == 1 and not node.keywords):
            return sp.sqrt(visit(node.args[0]))
        raise ValueError('Unsupported expression syntax in certificate')
    return visit(ast.parse(str(text), mode='eval'))


def geometric_expressions(payload, names):
    """Decode only the arithmetic entries of an exact coordinate certificate.

    A saved solver result is untrusted input. In particular, matrices must
    not be handed to SymPy while their entries are still expression strings.
    """
    if not isinstance(payload, dict):
        raise ValueError('A geometric certificate object is required')

    def matrix(value, label, rows=None):
        if (not isinstance(value, list) or (rows is not None and len(value) != rows)
                or any(not isinstance(row, list) or len(row) != 3 for row in value)):
            raise ValueError(f'Saved {label} must have three columns' +
                             (f' and {rows} rows' if rows is not None else ''))
        return [[expression(entry, names) for entry in row] for row in value]

    center = payload['center']
    if not isinstance(center, list) or len(center) != 3:
        raise ValueError('Saved center must be a three-vector')
    return {'metric': matrix(payload['metric'], 'metric', rows=3),
            'coordinates': matrix(payload['coordinates'], 'coordinates'),
            'center': [expression(entry, names) for entry in center]}


def verify_saved(folder):
    from verify import ExactField, verify_solution
    folder = Path(folder).resolve()
    case = json.loads((folder/'input.json').read_text())
    saved = json.loads((folder/'radii.json').read_text())
    alpha = sp.Symbol(saved['primitive_symbol'])
    names = {str(alpha): alpha}
    polynomial = expression(saved['primitive_polynomial'], names)
    interval = [expression(x, {}) for x in saved['isolating_interval']]
    radii = [expression(x, names) for x in saved['squared_radii']]
    if saved['groups'] != case['groups']:
        raise ValueError('Saved groups differ from the input partition')
    arguments = dict(primitive_symbol=alpha, primitive_polynomial=polynomial,
                     isolating_interval=interval)
    method = saved.get('verification_method', 'radial')
    if method == 'radial':
        report = verify_solution(case['faces'], case['groups'], radii, **arguments)
    elif method == 'exact_gram_geometry':
        from verify_gram import verify_gram_solution
        geometry = geometric_expressions(saved['geometric_certificate'], names)
        report = verify_gram_solution(case['faces'], case['groups'], radii,
                                      geometry, **arguments)
    else:
        raise ValueError('Unsupported saved verification method')
    if report.get('verified'):
        field = ExactField(alpha, polynomial, interval)
        T = sp.Symbol('T')
        expected = sp.Poly(sp.prod(T-x for x in radii), T).all_coeffs()
        actual = [expression(x, names) for x in saved['Q_coefficients']]
        if len(expected) != len(actual) or any(not field.zero(a-b) for a, b in zip(expected, actual)):
            report.update(verified=False, certified=False, status='UNVERIFIED',
                          reason='Saved Q coefficients do not equal the product of the radius factors')
    report['verification_method'] = method
    if report.get('verified') and (folder/'presentation.json').exists():
        from presentation_check import verify_presentation
        presentation = json.loads((folder/'presentation.json').read_text())
        presentation_report = verify_presentation(saved, presentation)
        report['presentation_check'] = presentation_report
        if presentation_report.get('verified') is not True:
            report.update(verified=False, certified=False, status='UNVERIFIED',
                          reason='The supplemental presentation does not match the exact certificate')
    return report


def worker(folder):
    try:
        report = verify_saved(folder)
    except Exception as exc:
        report = dict(verified=False, certified=False, status='UNVERIFIED',
                      reason=str(exc), exception=type(exc).__name__)
    # Coordinate verification includes exact plane covectors and offsets.
    # They are emitted as arithmetic text, like the solver's saved reports.
    print(json.dumps(report, indent=2, default=str))
    return 0 if report.get('verified') else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory')
    parser.add_argument('--time-limit', type=float, default=120.)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        return worker(args.output_directory)
    process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                                args.output_directory, '--worker'], start_new_session=True)
    try:
        return process.wait(timeout=args.time_limit)
    except (subprocess.TimeoutExpired, KeyboardInterrupt):
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
        print('UNVERIFIED: verification interrupted or exceeded its time limit')
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
