#!/usr/bin/env python3
"""BI numerical realization and conservative, interval-certified nonflatness screening."""

from proof import require
import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent


def signature(faces):
    def cycle(f):
        return min(tuple(f[i:] + f[:i]) for i in range(len(f)))

    return tuple(sorted(min(cycle(f), cycle(list(reversed(f)))) for f in faces))


def records(path):
    if path.suffix == ".jsonl":
        return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    value = json.loads(path.read_text())
    if isinstance(value, dict) and "faces" in value:
        return [value]
    if isinstance(value, list):
        if value and isinstance(value[0], list) and value[0] and isinstance(value[0][0], int):
            return [dict(faces=value)]
        return value
    raise ValueError("Expected face cycles, candidate objects, or JSONL")


def load(source, census=False):
    if census:
        cache = source
        for relative in ("validation/full-census", "cache_links_exact_extended"):
            if (source / relative).is_dir():
                cache = source / relative
                break
        files = sorted(cache.glob("candidates_*.jsonl"), key=lambda p: int(p.stem.split("_")[-1]))
        if not files:
            raise ValueError("No candidates_*.jsonl found")
        inputs = []
        for path in files:
            for value in records(path):
                item = dict(value)
                item["pass_number"] = int(path.stem.split("_")[-1])
                inputs.append(item)
    else:
        inputs = records(source)
    result = []
    seen = set()
    for rec in inputs:
        if not isinstance(rec, dict) or "faces" not in rec:
            raise ValueError("Candidate has no faces")
        # Only face cycles supply geometry. No accepted/rejected catalogs are read.
        from metric import validate_input

        validate_input(rec["faces"])
        sig = signature(rec["faces"])
        if sig in seen:
            continue
        seen.add(sig)
        item = {
            k: rec[k] for k in ("cid", "candidate_id", "name", "case", "pass_number") if k in rec
        }
        item["faces"] = rec["faces"]
        item.setdefault(
            "cid",
            item.get("candidate_id")
            or "metric-" + hashlib.sha256(repr(sig).encode()).hexdigest()[:16],
        )
        result.append(item)
    return result


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def decimals(interval, digits):
    from fractions import Fraction as Q

    lo, hi = map(Q, interval)
    scale = 10**digits
    a = lo.numerator * scale // lo.denominator
    b = -((-hi.numerator * scale) // hi.denominator)

    def render(n):
        sign = "-" if n < 0 else ""
        v = abs(n)
        return sign + str(v // scale) + "." + str(v % scale).zfill(digits)

    return "[" + render(a) + ", " + render(b) + "]"


def handoff(candidate, result, verified=None, numerical=None):
    record = dict(
        schema="crp-bi-numeric-handoff/2",
        cid=candidate["cid"],
        faces=candidate["faces"],
        status=result["status"],
        excluded=result.get("excluded", False),
        root_certified=result.get("root_certified", False),
        convex_metric_certified=result.get("convex_metric_certified", False),
        crp_membership_certified=False,
        exact_face_flatness_certified=False,
        needs_algebraic_stage=not result.get("excluded", False),
        reason=result.get("reason"),
    )
    if verified is not None:
        record.update(
            coordinate_intervals=verified["coordinate_intervals"],
            coordinate_frame=verified["coordinate_frame"],
            enclosure_target=verified["enclosure_target"],
            face_measurements=verified["faces"],
        )
    if numerical is not None and "positions" in numerical:
        record["bi_numerical_estimate"] = dict(
            positions=numerical["positions"],
            certified=False,
            frame="Raw BI frame, different from the gauged certificate frame",
            numerical_success=numerical.get("success", False),
            radii=numerical.get("radii"),
            squared_radii=numerical.get("squared_radii"),
            apex=numerical.get("apex"),
            vertex_orbits=numerical.get("vertex_orbits"),
        )
    return record


def attempt(candidate, directory, args):
    from metric import setup, initial_metric, orient_polygons, validate_input

    directory.mkdir()
    candidate = dict(candidate)
    faces = candidate["faces"]
    write_json(directory / "metadata.json", {k: v for k, v in candidate.items() if k != "faces"})
    result = dict(
        cid=candidate["cid"],
        status="UNRESOLVED_REALIZATION",
        excluded=False,
        root_certified=False,
        convex_metric_certified=False,
        crp_membership_certified=False,
        exact_face_flatness_certified=False,
        needs_algebraic_stage=True,
    )
    numerical = None
    verified = None
    start = time.monotonic()
    try:
        validate_input(faces)
        triangles, _ = initial_metric(faces)
        setup(dict(polygons=faces, flips=[], triangles=triangles, base=triangles[0]))
        faces = orient_polygons(faces)
        candidate["faces"] = faces
        write_json(directory / "input.json", dict(faces=faces))
        env = dict(
            os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1"
        )
        stages = [
            (
                "symmetry",
                [
                    sys.executable,
                    str(HERE / "symmetry.py"),
                    str(directory / "input.json"),
                    str(directory / "symmetry.json"),
                    str(directory / "metric.txt"),
                ],
            ),
            (
                "bi",
                [
                    str(HERE / "metric_probe"),
                    str(directory / "metric.txt"),
                    str(directory / "numerical.json"),
                    str(getattr(args, "max_iterations", 50000)),
                    str(getattr(args, "curvature_tolerance", 1e-9)),
                    str(getattr(args, "initial_radius_factor", 1.5)),
                    str(getattr(args, "radius_starts", 3)),
                ],
            ),
            (
                "proposal",
                [sys.executable, str(HERE / "produce.py"), str(directory), str(args.digits)],
            ),
        ]
        if getattr(args, "numerical_only", False):
            stages = stages[:2]
        for stage, command in stages:
            try:
                with (directory / (stage + ".log")).open("x") as log:
                    subprocess.run(
                        command,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                        timeout=args.timeout,
                        env=env,
                    )
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
                if stage == "proposal" and not (directory / "root-certificate.json").is_file():
                    # A timed-out proposal may have completed a Cartesian root
                    # inside one of its isolated routes. Replay it independently.
                    partial = sorted(directory.glob("proposal-*/root-certificate.json"))
                    if partial:
                        shutil.copyfile(partial[-1], directory / "root-certificate.json")
                if stage != "proposal" or not (directory / "root-certificate.json").is_file():
                    raise
                # Salvage only through a fresh independent replay, never by
                # trusting an incomplete producer report.
                (directory / "chosen-certificate.txt").write_text("root-certificate.json\n")
                result["producer_warning"] = (
                    "Extended proposal did not finish; independently replayed Cartesian certificate."
                )
            if stage == "bi":
                numerical = json.loads((directory / "numerical.json").read_text())[0]
                result["numerical_success"] = numerical.get("success", False)
                result["max_radial_defect"] = numerical.get("max_radial_defect")
                result["optimizer_reason"] = numerical.get("reason")
        if getattr(args, "numerical_only", False):
            result["status"] = (
                "NUMERICAL_ONLY" if numerical.get("success") else "UNRESOLVED_NUMERICAL"
            )
            result["reason"] = numerical.get("reason", "No numerical result")
            result["seconds"] = time.monotonic() - start
            write_json(directory / "result.json", result)
            write_json(directory / "handoff.json", handoff(candidate, result, None, numerical))
            return result
        chosen = (directory / "chosen-certificate.txt").read_text().strip()
        require(
            chosen in ("root-certificate.json", "hull-certificate.json"),
            "Unexpected certificate path",
        )
        cert_path = directory / chosen
        cert = json.loads(cert_path.read_text())
        require(cert["polygons"] == faces, "Certificate input differs from supplied metric")
        with (directory / "replay.log").open("x") as log:
            subprocess.run(
                [
                    sys.executable,
                    "-S",
                    str(HERE / "check.py"),
                    str(cert_path),
                    "--output",
                    str(directory / "verified.json"),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=args.timeout,
                env=env,
            )
        verified = json.loads((directory / "verified.json").read_text())
        require(verified["verified"] is True, "Certificate condition failed")
        for k in (
            "status",
            "excluded",
            "root_certified",
            "convex_metric_certified",
            "reason",
            "needs_algebraic_stage",
        ):
            result[k] = verified[k]
        result["certificate"] = cert_path.name
        lines = [
            "Coordinate interval target: " + verified["enclosure_target"],
            "Status: " + verified["status"],
            verified["reason"],
            "",
        ]
        for i, p in enumerate(verified["coordinate_intervals"]):
            lines.append(
                "v" + str(i) + " = (" + ", ".join(decimals(x, args.digits) for x in p) + ")"
            )
        (directory / "coordinates.txt").write_text("\n".join(lines) + "\n")
        write_json(directory / "face-measurements.json", verified["faces"])
        write_json(
            directory / "internal-folds.json",
            dict(note=verified["internal_fold_note"], measurements=verified["internal_folds"]),
        )
    except subprocess.TimeoutExpired as error:
        result["reason"] = "Stage timed out; no exclusion: " + Path(error.cmd[0]).name
    except subprocess.CalledProcessError as error:
        result["reason"] = "Stage failed; inspect logs. No exclusion: " + Path(error.cmd[0]).name
    except (ValueError, AssertionError, ArithmeticError, KeyError, IndexError, OSError) as error:
        result["reason"] = "Input or certification failure; no exclusion: " + str(error)
    result["seconds"] = time.monotonic() - start
    write_json(directory / "result.json", result)
    write_json(directory / "handoff.json", handoff(candidate, result, verified, numerical))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path)
    source.add_argument("--census", type=Path)
    select = ap.add_mutually_exclusive_group()
    select.add_argument("--id")
    select.add_argument("--case", type=int)
    select.add_argument("--all", action="store_true")
    select.add_argument("--list", action="store_true")
    ap.add_argument("--max-iterations", type=int, default=50000)
    ap.add_argument("--curvature-tolerance", type=float, default=1e-9)
    ap.add_argument("--initial-radius-factor", type=float, default=1.5)
    ap.add_argument(
        "--radius-starts",
        type=int,
        choices=(1, 2, 3),
        default=3,
        help="Bounded equal-radius starts: initial factor, 2.0, 1.2",
    )
    ap.add_argument(
        "--numerical-only", action="store_true", help="Skip certificate production and verification"
    )
    ap.add_argument("--digits", type=int, default=10)
    ap.add_argument("--timeout", type=float, default=300, help="Seconds per stage per candidate")
    ap.add_argument("--output-dir", type=Path)
    args = ap.parse_args()
    if not 3 <= args.digits <= 500:
        ap.error("--digits must be 3..500")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        ap.error("--timeout must be positive and finite")
    if args.max_iterations < 0:
        ap.error("--max-iterations must be nonnegative")
    if not math.isfinite(args.curvature_tolerance) or args.curvature_tolerance <= 0:
        ap.error("--curvature-tolerance must be positive and finite")
    if not math.isfinite(args.initial_radius_factor) or args.initial_radius_factor <= 1:
        ap.error("--initial-radius-factor must be finite and greater than one")
    candidates = load((args.census or args.input).expanduser().resolve(), bool(args.census))
    if args.list:
        for c in candidates:
            print(
                c["cid"],
                c.get("name", ""),
                f"V={len({v for f in c['faces'] for v in f})} F={len(c['faces'])}",
            )
        return 0
    if args.id:
        candidates = [
            c for c in candidates if args.id in (c.get("cid"), c.get("candidate_id"), c.get("name"))
        ]
    if args.case is not None:
        candidates = [c for c in candidates if c.get("case") == args.case]
    if not candidates:
        ap.error("No matching candidate")
    if not args.all and len(candidates) != 1:
        ap.error("Select --id, --case, or --all; inspect --list first")
    if not (HERE / "metric_probe").is_file():
        ap.error("Build the numerical helper with ./build.sh first")
    destination = (
        (args.output_dir or HERE / "results" / datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
        .expanduser()
        .resolve()
    )
    destination.mkdir(parents=True, exist_ok=False)
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob("*.py")}
    sources.update(
        {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for pattern in ("*.cpp", "*.h")
            for p in HERE.glob(pattern)
        }
    )
    sources["metric_probe"] = hashlib.sha256((HERE / "metric_probe").read_bytes()).hexdigest()
    write_json(
        destination / "manifest.json",
        dict(
            schema="crp-bi-certifier/4",
            source=str(args.census or args.input),
            count=len(candidates),
            digits=args.digits,
            timeout=args.timeout,
            sha256=sources,
            max_iterations=args.max_iterations,
            curvature_tolerance=args.curvature_tolerance,
            initial_radius_factor=args.initial_radius_factor,
            radius_starts=args.radius_starts,
            numerical_only=args.numerical_only,
            policy="Geometry comes only from regular-polygon face cycles. Full face-complex symmetry and orbit-constrained radial gradient descent. No algebraic coordinate recovery or CRP acceptance.",
        ),
    )
    reports = []
    (destination / "next-stage.jsonl").open("x").close()
    for i, c in enumerate(candidates):
        print(f"[{i+1}/{len(candidates)}] {c['cid']}", flush=True)
        slug = re.sub("[^a-zA-Z0-9_-]", "_", str(c["cid"]))[:100]
        report = attempt(c, destination / (str(i + 1).zfill(3) + "-" + slug), args)
        reports.append(report)
        if not report["excluded"]:
            case_dir = destination / (str(i + 1).zfill(3) + "-" + slug)
            with (destination / "next-stage.jsonl").open("a") as stream:
                stream.write(json.dumps(json.loads((case_dir / "handoff.json").read_text())) + "\n")
        print("  " + report["status"], flush=True)
        # Only the summary in this newly-created run is replaced as it progresses.
        (destination / "summary.json").write_text(json.dumps(reports, indent=2) + "\n")
    print("Results:", destination)
    return 1 if any(r["status"].startswith("UNRESOLVED") for r in reports) else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Interrupted; completed results are preserved.", file=sys.stderr)
        sys.exit(130)
