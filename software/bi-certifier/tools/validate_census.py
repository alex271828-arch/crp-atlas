#!/usr/bin/env python3
"""Run the new producer and independent checker on a saved candidate census.

This validation tool compares outcomes only after certification. Expected
outcomes never enter the solver. It does not rerun enumeration.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse complete validation records only with identical code and settings",
    )
    args = parser.parse_args()
    candidates = run.load(args.source, True)
    args.output.mkdir(parents=True, exist_ok=args.resume)
    settings = SimpleNamespace(
        digits=10,
        timeout=args.timeout,
        max_iterations=50000,
        curvature_tolerance=1e-9,
        initial_radius_factor=1.5,
        radius_starts=3,
        numerical_only=False,
    )
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("*.py")}
    hashes["metric_probe"] = hashlib.sha256((ROOT / "metric_probe").read_bytes()).hexdigest()
    manifest = dict(
        source=str(args.source),
        count=len(candidates),
        settings=vars(settings),
        sha256=hashes,
        workers=args.workers,
    )
    if args.resume:
        previous = json.loads((args.output / "manifest.json").read_text())
        for key in ("source", "count", "settings", "sha256"):
            if previous[key] != manifest[key]:
                raise ValueError("Cannot resume with changed " + key)
        run.write_json(args.output / f"resume-{time.time_ns()}.json", manifest)
    else:
        run.write_json(args.output / "manifest.json", manifest)

    def solve(item):
        i, candidate = item
        slug = re.sub("[^a-zA-Z0-9_-]", "_", str(candidate["cid"]))[:100]
        directory = args.output / f"{i:03d}-{slug}"
        if args.resume and directory.exists():
            result_path = directory / "result.json"
            if result_path.exists() and (directory / "handoff.json").exists():
                result = json.loads(result_path.read_text())
                supplied = json.loads((directory / "input.json").read_text())["faces"]
                if run.signature(supplied) != run.signature(candidate["faces"]):
                    raise ValueError("Saved validation input differs")
                if not result["status"].startswith("UNRESOLVED"):
                    result["directory"] = str(directory)
                    return result
            # Preserve interrupted/unresolved attempts outside the active case set.
            archive = args.output / "interrupted-attempts"
            archive.mkdir(exist_ok=True)
            directory.rename(archive / (directory.name + "-" + str(time.time_ns())))
        result = run.attempt(candidate, directory, settings)
        result["directory"] = str(directory)
        return result

    start = time.monotonic()
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(solve, enumerate(candidates, 1)):
            rows.append(row)
            print(len(rows), row["cid"], row["status"], round(row["seconds"], 3), flush=True)
            (args.output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    unchanged = all(
        hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == sha
        for name, sha in hashes.items()
    )
    run.write_json(
        args.output / "completion.json",
        dict(
            completed=True,
            count=len(rows),
            counts=dict(Counter(row["status"] for row in rows)),
            source_unchanged=unchanged,
            current_phase_wall_seconds=time.monotonic() - start,
            workers=args.workers,
            resumed=args.resume,
        ),
    )
    with (args.output / "next-stage.jsonl").open("x") as stream:
        for row in rows:
            if not row["excluded"]:
                stream.write(
                    json.dumps(json.loads((Path(row["directory"]) / "handoff.json").read_text()))
                    + "\n"
                )
    return int(not unchanged or any(row["status"].startswith("UNRESOLVED") for row in rows))


if __name__ == "__main__":
    sys.exit(main())
