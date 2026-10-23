#!/usr/bin/env python3
"""Replay saved certificates in isolated, standard-library-only processes."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=900)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    corpus = json.loads(args.corpus.read_text())

    def replay(item):
        source = Path(item["path"])
        name = source.parent.name
        output = args.output / (name + ".json")
        start = time.monotonic()
        with (args.output / (name + ".log")).open("x") as log:
            process = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    "-S",
                    str(args.module / "check.py"),
                    str(source),
                    "--output",
                    str(output),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
            )
        report = json.loads(output.read_text()) if output.exists() else {}
        return dict(
            case=name,
            certificate=str(source),
            sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            seconds=time.monotonic() - start,
            returncode=process.returncode,
            status=report.get("status"),
            expected=item["status"],
            match=process.returncode == 0 and report.get("status") == item["status"],
            convex_metric_certified=report.get("convex_metric_certified", False),
        )

    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = []
        for row in pool.map(replay, corpus):
            rows.append(row)
            print(len(rows), row["case"], row["status"], round(row["seconds"], 3), flush=True)
    summary = dict(
        total=len(rows),
        all_match=all(r["match"] for r in rows),
        counts=dict(Counter(r["status"] for r in rows)),
        wall_seconds=time.monotonic() - start,
        summed_case_seconds=sum(r["seconds"] for r in rows),
        workers=args.workers,
        rows=rows,
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0 if summary["all_match"] else 1


if __name__ == "__main__":
    sys.exit(main())
