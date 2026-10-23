#!/usr/bin/env python3
"""Select one of the 198 saved census survivors and invoke the v3 solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex

DEFAULT_PARTITIONS = None


def corrected_partition(path, data, corrections_directory=None):
    """Apply a packaged numerical-guide correction only to its exact source.

    Corrections change the input partition, never the algebraic equations or
    acceptance criteria. Hash checks bind them to the original partition and
    its BI coordinates and faces; custom partitions retain their own meaning.
    """
    directory = Path(corrections_directory or Path(__file__).resolve().parent/'partition_corrections')
    manifest_path = directory/'manifest.json'
    if not manifest_path.is_file():
        return path, data
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('format_version') != 1:
        raise ValueError('Unsupported partition-correction manifest version')
    source_hash = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    matches = [entry for entry in manifest['corrections']
               if entry['source_sha256'] == source_hash]
    if not matches:
        return path, data
    if len(matches) != 1:
        raise ValueError('More than one correction matches the selected partition')
    entry = matches[0]
    replacement = (directory/entry['partition']).resolve()
    if replacement.parent != directory.resolve():
        raise ValueError('Correction partition must be inside its catalog directory')
    corrected = json.loads(replacement.read_text())
    if any(corrected[key] != data[key] for key in ('census_index', 'display_id')):
        raise ValueError('Correction changes the selected body identity')
    bi = Path(data['input_file']).expanduser().resolve()
    if Path(corrected['input_file']).expanduser().resolve() != bi:
        raise ValueError('Correction changes the selected BI input')
    for source, expected in ((bi, entry['bi_sha256']),
                             (bi.parent/'input.json', entry['faces_sha256'])):
        if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise ValueError('Partition correction does not match the current BI coordinates or faces')
    if corrected.get('input_sha256') != entry['bi_sha256']:
        raise ValueError('Corrected partition has an inconsistent BI hash')
    if entry.get('replacement_sha256') and hashlib.sha256(replacement.read_bytes()).hexdigest() != entry['replacement_sha256']:
        raise ValueError('Packaged corrected partition has changed (SHA256 mismatch)')
    return replacement, corrected


def select_body(selector, directory):
    key = selector.strip().upper()
    matches = []
    for path in sorted(directory.glob('*.json')):
        data = json.loads(path.read_text())
        aliases = {str(data['census_index']), str(data['census_index']).zfill(3),
                   data['display_id'].upper(), 'CRP0-'+data['display_id'].upper(),
                   data['name'].upper()}
        if data.get('body_class') == 'Johnson':
            aliases.add(data['name'].split()[0].upper())
        if key in aliases:
            matches.append((path, data))
    if len(matches) != 1:
        raise ValueError(f'Expected one saved survivor matching {selector!r}; found {len(matches)}')
    return corrected_partition(*matches[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('body', help='CRP0-H identifier, original census number, or Johnson label such as J87')
    parser.add_argument('--partitions', type=Path, required=True, help='Directory of saved sphere partitions for this run')
    parser.add_argument('--out', type=Path, help='Output directory; default: package/out/<census number>-<CRP identifier>')
    parser.add_argument('--time-limit', type=float, default=120.)
    parser.add_argument('--radius-only', action='store_true', help='Skip the geometry solver; use only radius equations')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--presentation-time-limit', type=float, default=5.)
    parser.add_argument('--no-presentation', action='store_true')
    parser.add_argument('--show-command', action='store_true', help='Print the resolved command without running it')
    args = parser.parse_args()
    try:
        partition, data = select_body(args.body, args.partitions.expanduser().resolve())
    except (ValueError, KeyError, OSError) as exc:
        parser.error(str(exc))
    bi = Path(data['input_file']).expanduser().resolve()
    case = bi.parent
    for path in (bi, case/'input.json'):
        if not path.is_file():
            parser.error(f'Missing saved input: {path}')
    here = Path(__file__).resolve().parent
    body_id = 'CRP0-'+data['display_id']
    destination = args.out or here/'out'/f"{data['census_index']:03d}-{body_id}"
    command = [str(here/'run.sh'), '--case', str(case), '--bi', str(bi),
               '--spheres', str(partition), '--out', str(destination.expanduser().resolve()),
               '--time-limit', str(args.time_limit),
               '--presentation-time-limit', str(args.presentation_time_limit)]
    if args.overwrite:
        command.append('--overwrite')
    if args.radius_only:
        command.append('--radius-only')
    if args.no_presentation:
        command.append('--no-presentation')
    print(f"Census {data['census_index']}: {body_id} ({data['name']})", flush=True)
    if partition.parent.resolve() == (here/'partition_corrections').resolve():
        print(f"Using corrected sphere partition: {partition.name}", flush=True)
    if args.show_command:
        print(shlex.join(command))
        return 0
    os.execv(command[0], command)


if __name__ == '__main__':
    raise SystemExit(main())
