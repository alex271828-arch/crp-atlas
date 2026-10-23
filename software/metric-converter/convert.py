#!/usr/bin/env python3
"""Stage four: propose exact metric/Cartesian coordinates and replay separately."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from reconstruct import reconstruct, cartesian
from check import formula

ROOT = Path(__file__).resolve().parent


def write(path, data):
    with path.open("x") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def worker(source, destination):
    started = time.monotonic()
    case = json.loads((source / "input.json").read_text())
    saved = json.loads((source / "radii.json").read_text())
    A, G, points, center, meta = reconstruct(case, saved)
    B, cart, cc = cartesian(A, G, points, center)
    encode = lambda rows: [[A.encode(x) for x in row] for row in rows]
    source_keys = [
        "primitive_symbol",
        "primitive_polynomial",
        "isolating_interval",
        "groups",
        "squared_radii",
    ]
    # Exclude numerical positions, solver flags and pre-existing reports from proof data.
    cert = dict(
        schema="crp-metric-coordinates/1",
        input={k: case[k] for k in ["faces", "groups"]},
        source_algebra={k: saved[k] for k in source_keys},
        source_files={
            name: dict(
                path=str(source / name),
                sha256=hashlib.sha256((source / name).read_bytes()).hexdigest(),
            )
            for name in ["input.json", "radii.json"]
        },
        reconstruction=meta,
        metric=[[str(x) for x in row] for row in G.tolist()],
        positive_radicals=list(map(str, A.rad.generators)),
        metric_coordinates=encode(points),
        metric_center=[A.encode(x) for x in center],
        cartesian_transform=encode(B),
        cartesian_coordinates=encode(cart),
        cartesian_center=[A.encode(x) for x in cc],
    )
    for name, key in [
        ("metric_coordinate_formulas", "metric_coordinates"),
        ("cartesian_coordinate_formulas", "cartesian_coordinates"),
    ]:
        cert[name] = [[formula(x) for x in row] for row in cert[key]]
    write(destination / "coordinates.json", cert)
    write(
        destination / "proposal.json",
        dict(
            seconds=time.monotonic() - started, vertices=len(points), radicals=len(A.rad.generators)
        ),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--timeout", type=float, default=600)
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    if a.worker:
        worker(a.source.resolve(), a.output.resolve())
        return 0
    if not 0 < a.timeout < float("inf"):
        p.error("Timeout must be positive and finite")
    a.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    report = dict(status="UNVERIFIED", verified=False)
    try:
        with (a.output / "conversion.log").open("x") as log:
            run = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(ROOT / "convert.py"),
                    str(a.source.resolve()),
                    "--output",
                    str(a.output.resolve()),
                    "--worker",
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=a.timeout,
            )
        if run.returncode:
            raise RuntimeError("Coordinate proposal failed; see conversion.log")
        with (a.output / "verification.log").open("x") as log:
            run = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(ROOT / "check.py"),
                    str((a.output / "coordinates.json").resolve()),
                    "--output",
                    str((a.output / "verification.json").resolve()),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=a.timeout,
            )
        report = json.loads((a.output / "verification.json").read_text())
        if run.returncode:
            report["verified"] = False
            report["status"] = "UNVERIFIED"
    except (subprocess.TimeoutExpired, RuntimeError) as exc:
        report["reason"] = str(exc)

    if report["verified"]:
        from report import render

        render(a.output)
    report["seconds"] = time.monotonic() - start
    write(a.output / "result.json", report)
    print(json.dumps(report, indent=2))
    return 0 if report["verified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
