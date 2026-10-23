#!/usr/bin/env python3
"""Run the ten worked examples with independent wall-clock limits."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', default='out/worked-examples')
    parser.add_argument('--time-limit', type=float, default=120.)
    parser.add_argument('--radius-only', action='store_true', help='Skip the geometry solver; use only radius equations')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--presentation-time-limit', type=float, default=5.)
    parser.add_argument('--no-presentation', action='store_true')
    parser.add_argument('--recheck', action='store_true', help='Also run the standalone saved-certificate checker')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    started = time.monotonic()
    for name in ('CRP87', 'CRP140', 'J20', 'J84', 'J85', 'J86', 'J87', 'J88', 'J89', 'J92'):
        case, dest = here/'examples'/name, out/name
        command = [sys.executable, str(here/'solve.py'), '--case', str(case),
                   '--spheres', str(case/'spheres.json'), '--out', str(dest),
                   '--time-limit', str(args.time_limit),
                   '--presentation-time-limit', str(args.presentation_time_limit)]
        if args.overwrite:
            command.append('--overwrite')
        if args.radius_only:
            command.append('--radius-only')
        if args.no_presentation:
            command.append('--no-presentation')
        run = subprocess.run(command)
        result_path = dest/'result.json'
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        row = {'case': name, 'status': result.get('status', 'NO_RESULT'),
               'seconds': result.get('seconds'), 'exit_code': run.returncode,
               'primitive_polynomial': result.get('primitive_polynomial'),
               'Q_coefficients': result.get('Q_coefficients'),
               'presentation_status': result.get('presentation_status'),
               'presentation_seconds': result.get('presentation_seconds'),
               'presentation_improved': result.get('presentation_improved'),
               'output_directory': str(dest)}
        if args.recheck and run.returncode == 0 and result.get('status') == 'SOLVED':
            with (dest/'independent-check.json').open('w') as log:
                check = subprocess.run([sys.executable, str(here/'check.py'), str(dest),
                                       '--time-limit', str(args.time_limit)], stdout=log)
            row['independent_check_passed'] = check.returncode == 0
        rows.append(row)
        print(f"{name}: {row['status']} ({row['seconds']} seconds)", flush=True)
    report = {'cases': rows, 'total_seconds': time.monotonic()-started,
              'all_solved': all(r['status'] == 'SOLVED' and r['exit_code'] == 0 for r in rows),
              'all_presentations_verified': (None if args.no_presentation else
                                             all(r.get('presentation_status') == 'VERIFIED' for r in rows)),
              'all_independent_checks_passed': all(r.get('independent_check_passed', not args.recheck) for r in rows)}
    (out/'benchmark.json').write_text(json.dumps(report, indent=2)+'\n')
    return 0 if report['all_solved'] and report['all_independent_checks_passed'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
