#!/usr/bin/env python3
"""Produce an optional compact presentation from an existing exact certificate."""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


def save(path, data):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, default=str, allow_nan=False)+'\n')
    temporary.replace(path)


def baseline(saved):
    """An identity presentation remains available even if simplification times out."""
    symbol = safe_symbol_name(saved['primitive_symbol'])
    polynomial = str(saved['primitive_polynomial'])
    values = [str(x) for x in saved['squared_radii']]
    size = len(polynomial)+sum(map(len, values))
    return dict(format_version=1, verified=True, method='original_certificate',
                generator=dict(symbol=symbol, polynomial=polynomial,
                    isolating_interval=[str(x) for x in saved['isolating_interval']],
                    definition_in_source=symbol,
                    description='Original certified primitive element'),
                definitions=[], squared_radii=values,
                Q_factored='*'.join(f'(T-s{i})' for i in range(len(values))),
                score=size, baseline_score=size, seconds=0.)


def safe_symbol_name(value):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', value):
        raise ValueError('A report symbol must be an ASCII arithmetic identifier')
    if value == 'T':
        raise ValueError('T is reserved for the radius polynomial')
    return value


def mathml_fragment(value):
    """Emit well-formed MathML Core, with explicit fences for current browsers."""
    import sympy as sp
    text = sp.printing.mathml(value, printer='presentation')
    # SymPy emits literal comparison signs and HTML-only named entities.
    text = text.replace('<mo><</mo>', '<mo>&lt;</mo>').replace('<mo>></mo>', '<mo>&gt;</mo>')
    text = re.sub(r'&([A-Za-z][A-Za-z0-9]+);',
                  lambda m: m.group(0) if m[1] in ('lt', 'gt', 'amp', 'quot', 'apos')
                  else html.unescape(m.group(0)), text)
    root = ET.fromstring(text)
    for node in root.iter():
        if node.tag != 'mfenced':
            continue
        opening, closing = node.get('open', '('), node.get('close', ')')
        separators = ''.join(node.get('separators', ',').split())
        children = list(node)
        node.clear()
        node.tag = 'mrow'
        if opening:
            ET.SubElement(node, 'mo', fence='true').text = opening
        for i, child in enumerate(children):
            if i and separators:
                ET.SubElement(node, 'mo', separator='true').text = separators[min(i-1, len(separators)-1)]
            node.append(child)
        if closing:
            ET.SubElement(node, 'mo', fence='true').text = closing
    return ET.tostring(root, encoding='unicode')


def render(presentation, saved, out, name):
    """Write readable Markdown and offline MathML, using the same expressions."""
    import sympy as sp
    from check import expression
    out = Path(out)
    generator = presentation['generator']
    t = sp.Symbol(safe_symbol_name(generator['symbol']))
    names = {str(t): t}
    polynomial = expression(generator['polynomial'], names)
    interval = [expression(x, {}) for x in generator['isolating_interval']]
    rows = [('Defining equation', sp.Eq(polynomial, 0, evaluate=False))]
    if interval[0] == interval[1]:
        root = sp.Eq(t, interval[0], evaluate=False)
    else:
        root = sp.And(sp.Lt(interval[0], t, evaluate=False),
                      sp.Lt(t, interval[1], evaluate=False), evaluate=False)
    rows.append(('Chosen real root', root))
    for definition in presentation['definitions']:
        symbol = sp.Symbol(safe_symbol_name(definition['symbol']))
        value = expression(definition['expression'], names)
        rows.append(('Shared quantity', sp.Eq(symbol, value, evaluate=False)))
        names[str(symbol)] = symbol
    for i, text in enumerate(presentation['squared_radii']):
        value = expression(text, names)
        rows.append((f'Squared radius {i}', sp.Eq(sp.Symbol(f's_{i}'), value, evaluate=False)))
    T = sp.Symbol('T')
    factors = sp.Mul(*(T-sp.Symbol(f's_{i}') for i in range(len(presentation['squared_radii']))),
                     evaluate=False)
    rows.append(('Radius polynomial', sp.Eq(sp.Function('Q')(T), factors, evaluate=False)))
    description = generator.get('description', '')
    improved = presentation['score'] < presentation['baseline_score']
    explanation = ('A shorter presentation was found.' if improved else
                   'The original presentation was retained; no shorter verified form was found within the budget.')
    metrics = (f"Formula size: {presentation['score']} characters; "
               f"baseline: {presentation['baseline_score']}. "
               'The score counts the defining polynomial, shared definitions, and squared-radius formulas; '
               'certificate conversion maps are retained separately.')
    md = [f'# {name}: exact radii', '', explanation, '', description, '', metrics, '']
    body = []
    for label, value in rows:
        md.extend([f'**{label}**', '', '\\[', sp.latex(value), '\\]', ''])
        mathml = mathml_fragment(value)
        if label == 'Radius polynomial':
            # SymPy lowercases an undefined function's MathML name.
            mathml = '<mrow><mi>Q</mi><mo>(</mo><mi>T</mi><mo>)</mo><mo>=</mo>'+mathml_fragment(factors)+'</mrow>'
        body.append(f'<section><h2>{html.escape(label)}</h2><div class="formula">'
                    f'<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">{mathml}</math>'
                    '</div></section>')
    md.extend(['The radii are the positive square roots: \\(r_i=\\sqrt{s_i}\\).', '',
               'Vertex indices use the original input numbering.', ''])
    groups = []
    for i, vertices in enumerate(saved['groups']):
        md.append(f'- Sphere {i}: vertices '+', '.join(map(str, vertices)))
        groups.append(f'<li>Sphere {i}: '+html.escape(', '.join(map(str, vertices)))+'</li>')
    md.extend(['', 'Every displayed formula is an exact expression for the corresponding certified squared radius.',
               'The original field representation and geometric certificate remain in `radii.json` and `certificate.json`.',
               'The exact change of generator is recorded in `presentation.json`.', ''])
    page = '''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Exact radii</title><style>
body{font:17px/1.6 system-ui,sans-serif;color:#172330;background:#f7f9fb;margin:0}
main{max-width:1100px;margin:40px auto;padding:0 28px 40px}h1{font-size:30px}
h2{font-size:16px;color:#486278;margin:0 0 14px}section{background:white;
border:1px solid #dce4ec;border-radius:10px;padding:20px;margin:18px 0}
.formula{overflow-x:auto;padding:10px 0}math{font-size:19px}p{max-width:900px}
.muted{color:#586b7e;font-size:14px}a{color:#175e9e}li{overflow-wrap:anywhere}
</style><main>'''
    page += f'<h1>{html.escape(name)}: exact radii</h1><p>{html.escape(explanation)}</p>'
    page += f'<p>{html.escape(description)}</p><p class="muted">{html.escape(metrics)}</p>'
    page += ''.join(body)
    page += '<p>The radii are the positive square roots of the squared radii above.</p>'
    page += '<h2>Vertex membership</h2><p class="muted">Original input numbering.</p><ul>'+''.join(groups)+'</ul>'
    page += ('<p class="muted">Every displayed formula equals the corresponding certified squared radius exactly. '
             'The full change of generator is in <a href="presentation.json">presentation.json</a>.</p></main></html>')
    (out/'presentation.md').write_text('\n'.join(md))
    (out/'presentation.html').write_text(page)


def worker(source, destination, seconds):
    from presentation import build_presentation
    from presentation_check import verify_presentation
    source, destination = Path(source), Path(destination)
    saved = json.loads((source/'radii.json').read_text())
    result = json.loads((source/'result.json').read_text())
    provenance = {'radii_sha256': hashlib.sha256((source/'radii.json').read_bytes()).hexdigest(),
                  'source_directory': str(source), 'solver_version': result.get('version')}
    derivation_file = source/'derivation.json'
    derivation = json.loads(derivation_file.read_text()) if derivation_file.exists() else None
    last_signature, last_report = None, None

    def publish(candidate):
        nonlocal last_signature, last_report
        signature = json.dumps({key: candidate[key] for key in
            ('format_version', 'generator', 'definitions', 'squared_radii', 'Q_factored')}, sort_keys=True)
        report = (last_report if signature == last_signature else
                  verify_presentation(saved, candidate))
        if report.get('verified') is not True:
            raise ValueError('Presentation failed independent equality check: '+str(report))
        last_signature, last_report = signature, report
        candidate = dict(candidate, source=provenance, verification=report)
        save(destination/'presentation.json', candidate)
        print(json.dumps({'method': candidate['method'], 'score': candidate['score'],
                          'baseline_score': candidate['baseline_score']}), flush=True)

    candidate = build_presentation(saved, derivation=derivation,
                                   max_seconds=max(.01, seconds-.5), on_improvement=publish)
    publish(candidate)
    return 0


def run_presentation(source, destination=None, *, time_limit=5.):
    source = Path(source).resolve()
    destination = Path(destination or source).resolve()
    result = json.loads((source/'result.json').read_text())
    if result.get('status') != 'SOLVED' or result.get('certified') is not True:
        raise ValueError('A SOLVED, certified source result is required')
    saved = json.loads((source/'radii.json').read_text())
    destination.mkdir(parents=True, exist_ok=True)
    initial = baseline(saved)
    initial['source'] = {'radii_sha256': hashlib.sha256((source/'radii.json').read_bytes()).hexdigest(),
                         'source_directory': str(source), 'solver_version': result.get('version')}
    save(destination/'presentation.json', initial)
    start = time.monotonic()
    timed_out = False
    with (destination/'presentation.log').open('w') as log:
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), str(source),
                                    '--out', str(destination), '--worker', '--time-limit', str(time_limit)],
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            process.wait(timeout=time_limit)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            if isinstance(exc, KeyboardInterrupt):
                raise
            timed_out = True
    elapsed = time.monotonic()-start
    candidate = json.loads((destination/'presentation.json').read_text())
    candidate.update(search_seconds=elapsed,
                     time_limit_reached=timed_out or candidate.get('budget_exhausted', False),
                     hard_timeout=timed_out,
                     search_exit_code=process.returncode)
    save(destination/'presentation.json', candidate)
    render(candidate, saved, destination, result.get('name', 'Polyhedron'))
    return {'presentation_status': 'VERIFIED',
            'presentation_improved': candidate['score'] < candidate['baseline_score'],
            'presentation_seconds': elapsed,
            'presentation_time_limit_reached': candidate['time_limit_reached'],
            'presentation_files': ['presentation.html', 'presentation.md', 'presentation.json']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', help='Directory containing an existing certified result')
    parser.add_argument('--out', help='Presentation destination; defaults to the source directory')
    parser.add_argument('--time-limit', type=float, default=5., help='Separate simplification budget in seconds')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not math.isfinite(args.time_limit) or args.time_limit <= 0:
        parser.error('Presentation time limit must be finite and positive')
    if args.worker:
        return worker(args.source, args.out or args.source, args.time_limit)
    report = run_presentation(args.source, args.out, time_limit=args.time_limit)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
