#!/usr/bin/env python3
"""Compare every candidate complex against a previous run; no catalogue is used."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from link_growth import canonical_code
from storage import rows, atomic_json


def compare(reference, current):
    report = []
    for number in range(63):
        old_file = reference / f"pass_{number}_done.json"
        new_file = current / f"pass_{number}_done.json"
        before = json.loads(old_file.read_text()) if old_file.exists() else {}
        after = json.loads(new_file.read_text()) if new_file.exists() else {}
        old_rows = rows(reference / f"candidates_{number}.jsonl")
        new_rows = rows(current / f"candidates_{number}.jsonl")
        old_codes = {canonical_code(r["faces"]) for r in old_rows}
        new_codes = {canonical_code(r["faces"]) for r in new_rows}

        def identifiers(codes):
            return sorted(hashlib.sha256(repr(c).encode()).hexdigest() for c in codes)

        missing = identifiers(old_codes - new_codes)
        extra = identifiers(new_codes - old_codes)
        complete = before.get("enumeration_exhausted", False) and after.get(
            "enumeration_exhausted", False
        )
        report.append(
            dict(
                pass_number=number,
                complete=complete,
                previous=len(old_rows),
                current=len(new_rows),
                missing=missing,
                extra=extra,
                duplicates=len(new_rows) - len(new_codes),
                match=complete
                and not missing
                and not extra
                and len(new_rows) == len(new_codes),
            )
        )
    return dict(
        all_passes_match=all(r["match"] for r in report),
        completed_passes=sum(r["complete"] for r in report),
        previous_candidates=sum(r["previous"] for r in report),
        current_candidates=sum(r["current"] for r in report),
        passes=report,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("current", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.reference, args.current)
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != "passes"}, indent=2))
    raise SystemExit(0 if result["all_passes_match"] else 1)
