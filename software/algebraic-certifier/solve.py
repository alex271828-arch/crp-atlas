#!/usr/bin/env python3
"""Exact squared-radius solver. Run --help for the public interface."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

VERSION = '3.0.11'


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [jsonable(v) for v in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def save(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(jsonable(value), indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def progress(out, stage, **details):
    save(out/'progress.json', dict(stage=stage, **details))
    print(stage, flush=True)


def geometry_search_budget(remaining_seconds, radial_applicable):
    """Reserve radial fallback time only when that route can be used.

    The caller's remaining budget already excludes the final verification and
    output reserve; the parent process still enforces the total wall-clock limit.
    """
    if not radial_applicable:
        return remaining_seconds
    return min(65., remaining_seconds*.65)


def certify_candidates(case, solution, out, base, reports, route):
    """Accept a search route only after checking the original geometry."""
    import sympy as sp
    from verify import verify_solution, ExactField
    candidates = solution.get('candidates', [])
    if not candidates and solution.get('radius_expressions'):
        candidates = [solution]
    for index, candidate in enumerate(candidates):
        progress(out, 'verifying_candidate', route=route,
                 candidate=index+1, candidates=len(candidates))
        expressions = candidate.get('squared_radii', candidate.get('radius_expressions'))
        if isinstance(expressions, dict):
            symbols = sp.symbols(f"s0:{len(case['groups'])}", positive=True)
            expressions = [expressions[v] if v in expressions else expressions[str(v)] for v in symbols]
        field_arguments = dict(primitive_symbol=candidate.get('primitive_symbol'),
                    primitive_polynomial=candidate.get('primitive_polynomial'),
                    isolating_interval=candidate.get('isolating_interval'))
        geometry = candidate.get('geometric_certificate')
        verification_method = 'radial'
        radial_applicable = case.get('preflight', {}).get('center_clearance', 1.) > 1e-8
        if radial_applicable:
            report = verify_solution(case['faces'], case['groups'], expressions, **field_arguments)
        else:
            report = dict(verified=False, certified=False, status='UNVERIFIED',
                          reason='An arbitrary-center sphere cover needs an exact geometric certificate')
        if not report.get('verified') and geometry is not None:
            from verify_gram import verify_gram_solution
            radial_report = report
            report = verify_gram_solution(case['faces'], case['groups'], expressions,
                                         geometry, **field_arguments)
            report['radial_attempt'] = radial_report
            verification_method = 'exact_gram_geometry'
        report['search_route'] = route
        reports.append(report)
        save(out/'verification.json', reports)
        if report.get('verified') is not True:
            continue
        T = sp.Symbol('T')
        polynomial = sp.Poly(sp.prod(T-s for s in expressions), T)
        field = ExactField(candidate.get('primitive_symbol'),
                           candidate.get('primitive_polynomial'),
                           candidate.get('isolating_interval'))
        coefficients = [field.normal(c) for c in polynomial.all_coeffs()]
        radii = {'primitive_symbol': candidate.get('primitive_symbol'),
                 'primitive_polynomial': candidate.get('primitive_polynomial'),
                 'isolating_interval': candidate.get('isolating_interval'),
                 'groups': case['groups'], 'squared_radii': expressions,
                 'radii': [sp.sqrt(s) for s in expressions],
                 'Q_factored': '*'.join(f'(T-({s}))' for s in expressions),
                 'Q_coefficients': coefficients,
                 'verification_method': verification_method,
                 'auxiliary_assignments': candidate.get('auxiliary_assignments', {})}
        # Preserve the candidate's original expression structure as a display
        # hint. Canonical field elements remain the certificate authority.
        radii['presentation_hints'] = {
            'squared_radii': expressions,
            'auxiliary_assignments': candidate.get('auxiliary_assignments', {}),
            'search_route': route}
        if verification_method == 'exact_gram_geometry':
            radii['geometric_certificate'] = geometry
            base['center_location'] = report['center_location']
        save(out/'radii.json', radii)
        save(out/'certificate.json', report)
        base.update(status='SOLVED', certified=True, search_route=route,
                    verification_method=verification_method,
                    squared_radii=expressions, Q_coefficients=coefficients,
                    primitive_polynomial=candidate.get('primitive_polynomial'),
                    isolating_interval=candidate.get('isolating_interval'),
                    reason=('Exact coordinates, regular faces, convexity, and sphere membership verified'
                            if verification_method == 'exact_gram_geometry' else
                            'Exact algebraic solution verified against radial closure, flatness, and convexity'))
        return True
    return False


def worker(args):
    started = time.monotonic()
    out = Path(args.out).resolve()
    radius_only = getattr(args, 'radius_only', False)
    base = {'version': VERSION, 'status': 'UNRESOLVED', 'certified': False,
            'radius_only': radius_only}
    reports = []

    def finish():
        save(out/'result.json', dict(base, seconds=time.monotonic()-started))
        return 0 if base['status'] in ('SOLVED', 'EQUATIONS_ONLY') else 2

    def remaining():
        return max(.1, args.time_limit-(time.monotonic()-started)-2.)

    def accept(solution, route):
        save(out/'derivation.json', solution)
        return certify_candidates(case, solution, out, base, reports, route)

    try:
        from inputs import load_case
        progress(out, 'reading_input')
        case = load_case(args.case, bi=args.bi, spheres=args.spheres, tolerance=args.tolerance)
        save(out/'input.json', case)
        base.update(name=case['name'], sources=case['sources'], notes=case['notes'])
        radial_applicable = case['preflight']['center_clearance'] > 1e-8
        if not radial_applicable:
            if radius_only:
                base.update(reason='Radius-only mode requires a strictly interior sphere center; the supplied center does not pass the numerical interior check',
                            algebraic_status='UNRESOLVED', candidates_checked=0)
                return finish()
            case['notes'].append('Using exact metric geometry for a sphere center that is not numerically strictly interior')
            save(out/'input.json', case)
        if case['preflight']['face_pair_error'] > 1e-4:
            raise ValueError('BI coordinates are too far from regular flat faces to guide this solver (face distance error > 1e-4)')
        from encoder import build_system
        from patches import enrich_system
        from reduce import solve_model, _time_limit, ReductionTimeout

        def radius_model():
            result = build_system(case['faces'], case['groups'], case['squared_radii_numeric'],
                                  defer_dense=True)
            enrich_system(result)
            return result

        if args.equations_only:
            if radial_applicable:
                progress(out, 'generating_exact_radius_equations')
                save(out/'equations.json', radius_model())
            else:
                from gram_model import build_gram_model, symmetry_blocks
                from distance_kernel import gram_rank_relations
                progress(out, 'generating_exact_metric_equations')
                metric_model = build_gram_model(case)
                blocks = symmetry_blocks(metric_model, positions=case['positions'])
                relations = [gram_rank_relations(b['gram'], b['numerical_rank'],
                             witness=metric_model['witness']) for b in blocks]
                save(out/'equations.json', dict(method='exact_metric',
                    gram=metric_model['gram'].tolist(), variables=metric_model['variables'],
                    block_rank_relations=relations, derivation=metric_model['trace']))
            base.update(status='EQUATIONS_ONLY', reason='Exact equations generated; solving was not requested')
            return finish()
        model = None
        solution = {'status': 'UNRESOLVED', 'candidates': []}
        if radial_applicable:
            progress(out, 'trying_small_local_radius_system')
            try:
                with _time_limit(min(3., max(.1, args.time_limit*.1))):
                    model = radius_model()
                    solution = solve_model(model, max_seconds=min(2., args.time_limit*.07), cheap_only=True)
            except ReductionTimeout:
                model = None
                solution = {'status': 'UNRESOLVED', 'candidates': [],
                            'reason': 'Small local radius prepass reached its budget'}
        else:
            solution = {'status': 'UNRESOLVED', 'candidates': [],
                        'reason': 'Radial prepass skipped: the supplied sphere center is not strictly interior'}
        save(out/'local-derivation.json', solution)
        if solution.get('candidates'):
            save(out/'equations.json', model)
        if accept(solution, 'local_radius'):
            return finish()

        if not radius_only:
            from gram_solve import solve_gram_model
            geometry_seconds = geometry_search_budget(remaining(), radial_applicable)
            progress(out, 'solving_exact_distance_and_symmetry_equations',
                     search_budget_seconds=geometry_seconds,
                     radial_fallback_available=radial_applicable)
            solution = solve_gram_model(case, max_seconds=geometry_seconds)
            save(out/'gram-derivation.json', solution)
            save(out/'equations.json', {'method': 'exact_invariant_gram_completion',
                                       'derivation': solution.get('derivation', [])})
            if accept(solution, 'symmetry_gram'):
                return finish()

        if not radial_applicable:
            base.update(reason=('No exact metric candidate passed geometric verification' if reports else
                                solution.get('reason', 'No isolated exact metric candidate found')),
                        algebraic_status=solution.get('status'), candidates_checked=len(reports))
            return finish()

        progress(out, 'generating_local_radius_relations')
        if model is None:
            model = radius_model()
        save(out/'equations.json', model)
        progress(out, 'solving_radius_and_short_angle_equations')
        solution = solve_model(model, max_seconds=remaining()*.8)
        save(out/'radial-derivation.json', solution)
        if accept(solution, 'radial_only' if radius_only else 'radial_fallback'):
            return finish()
        base.update(reason=('No candidate passed exact geometric verification' if reports else
                            solution.get('reason', 'No isolated algebraic candidate found')),
                    algebraic_status=solution.get('status'), candidates_checked=len(reports))
        return finish()
    except Exception as exc:
        traceback.print_exc()
        base.update(reason=str(exc), exception=type(exc).__name__)
        return finish()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=f'%(prog)s {VERSION}')
    parser.add_argument('--case', required=True, help='BI case directory, or census JSON containing faces')
    parser.add_argument('--bi', help='BI numerical.json; defaults to the case directory')
    parser.add_argument('--spheres', help='Existing numerical sphere partition JSON; otherwise discover a partition')
    parser.add_argument('--out', required=True, help='Directory for equations, derivation, result, and certificate')
    parser.add_argument('--time-limit', type=float, default=120., help='Hard wall-clock limit in seconds (default: 120)')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='Numerical sphere clustering tolerance only')
    parser.add_argument('--equations-only', action='store_true', help='Save exact equations without solving')
    parser.add_argument('--radius-only', action='store_true',
                        help='Use only radius equations, skipping the geometry solver (requires a strictly interior sphere center)')
    parser.add_argument('--overwrite', action='store_true', help='Allow replacing a previous run in the output directory')
    parser.add_argument('--presentation-time-limit', type=float, default=5.,
                        help='Separate optional simplification budget after certification (default: 5 seconds)')
    parser.add_argument('--no-presentation', action='store_true',
                        help='Write only the original exact certificate files')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if (not math.isfinite(args.time_limit) or not math.isfinite(args.tolerance)
            or not math.isfinite(args.presentation_time_limit)
            or args.time_limit <= 0 or args.tolerance <= 0 or args.presentation_time_limit <= 0):
        parser.error('Time limits and tolerance must be finite and positive')
    out = Path(args.out).expanduser().resolve()
    args.out = str(out)
    out.mkdir(parents=True, exist_ok=True)
    if args.worker:
        return worker(args)
    if (out/'result.json').exists() and not args.overwrite:
        parser.error('Output already contains a result; choose another directory or use --overwrite')
    if args.overwrite:
        for name in ('result.json', 'radii.json', 'certificate.json', 'verification.json',
                     'derivation.json', 'equations.json', 'progress.json', 'input.json',
                     'independent-check.json', 'gram-derivation.json', 'local-derivation.json', 'radial-derivation.json',
                     'presentation.json', 'presentation.md', 'presentation.html', 'presentation.log', 'presentation.json.tmp'):
            (out/name).unlink(missing_ok=True)
    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--case', args.case,
               '--out', str(out), '--time-limit', str(args.time_limit), '--tolerance', str(args.tolerance)]
    if args.bi:
        command += ['--bi', args.bi]
    if args.spheres:
        command += ['--spheres', args.spheres]
    if args.equations_only:
        command += ['--equations-only']
    if args.radius_only:
        command += ['--radius-only']
    start = time.monotonic()
    try:
        with (out/'run.log').open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            try:
                process.wait(timeout=args.time_limit)
            except subprocess.TimeoutExpired:
                import signal
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                stage = json.loads((out/'progress.json').read_text()) if (out/'progress.json').exists() else {}
                save(out/'result.json', {'version': VERSION, 'status': 'UNRESOLVED', 'certified': False,
                     'radius_only': args.radius_only,
                     'reason': 'Hard wall-clock limit reached', 'last_progress': stage,
                     'seconds': time.monotonic()-start})
    except KeyboardInterrupt:
        if 'process' in locals() and process.poll() is None:
            import signal
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        save(out/'result.json', {'version': VERSION, 'status': 'UNRESOLVED', 'certified': False,
                                'radius_only': args.radius_only,
                                'reason': 'Interrupted by user', 'seconds': time.monotonic()-start})
        return 130
    if not (out/'result.json').exists():
        save(out/'result.json', {'version': VERSION, 'status': 'UNRESOLVED', 'certified': False,
                                'radius_only': args.radius_only,
                                'reason': f'Worker exited with code {process.returncode}; inspect run.log'})
    result = json.loads((out/'result.json').read_text())
    if result.get('status') == 'SOLVED' and result.get('certified') and not args.no_presentation:
        try:
            from present import run_presentation
            result.update(run_presentation(out, time_limit=args.presentation_time_limit))
        except KeyboardInterrupt:
            result.update(presentation_status='INTERRUPTED')
            save(out/'result.json', result)
            print('Presentation interrupted; the completed exact certificate is saved.')
            return 130
        except Exception as exc:
            # A display limitation must not change a completed proof result.
            result.update(presentation_status='UNAVAILABLE', presentation_reason=str(exc))
        save(out/'result.json', result)
    print(f"{result['status']}: {result.get('reason', '')}")
    print(f"Output: {out}")
    if result.get('presentation_status') == 'VERIFIED':
        print(f"Readable radii: {out/'presentation.html'}")
    return 0 if result['status'] in ('SOLVED', 'EQUATIONS_ONLY') else 2


if __name__ == '__main__':
    raise SystemExit(main())
