#!/usr/bin/env python3
"""Census v1.16: enumerate CRP° candidates for the BI and algebraic stages."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import signal
import subprocess
import sys
import time

import census_alphabet as alphabet
from face_types import FACE_SIZES, ANGLE_UNIT_DEGREES
import link_engine as engine
import link_growth as growth
from storage import append, atomic_json, lock, rows

HERE = Path(__file__).resolve().parent
VERSION = "1.16"


def cache(root=None):
    directory = Path(root).resolve() if root else HERE / "cache"
    directory.mkdir(parents=True, exist_ok=True)
    sources = sorted(p for p in HERE.glob("*.py") if not p.name.startswith("test_"))
    configuration = dict(
        version=VERSION,
        face_sizes=list(FACE_SIZES),
        angle_unit_degrees=ANGLE_UNIT_DEGREES,
        algorithm="all-link-completions",
        coplanar_faces=True,
        hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    )
    path = directory / "configuration.json"
    with lock(directory / "configuration.lock", blocking=True):
        if path.exists():
            if json.loads(path.read_text()) != configuration:
                raise ValueError(
                    "Cache source/configuration mismatch; use a fresh --cache-dir"
                )
        else:
            atomic_json(path, configuration)
    return directory


def summary(directory, pass_number, state, candidate_count):
    record = dict(
        version=VERSION,
        face_sizes=FACE_SIZES,
        pass_number=pass_number,
        link=alphabet.ORDER[pass_number],
        enumeration_exhausted=state["exhausted"],
        geometric_certification_complete=False,
        candidates=candidate_count,
        stats=state["stats"],
        explanation="All candidates require the separate BI and algebraic stages.",
    )
    atomic_json(directory / f"pass_{pass_number}_summary.json", record)
    if state["exhausted"]:
        atomic_json(directory / f"pass_{pass_number}_done.json", record)
    return record


def grow(directory, pass_number, seconds, stopped=lambda: False):
    if not 0 <= pass_number < len(alphabet.SIGS):
        raise ValueError("Pass number is out of range")
    with lock(directory / f"pass_{pass_number}.lock"):
        checkpoint = directory / f"pass_{pass_number}_state.pkl"
        menu = tuple(alphabet.SIGS[pass_number:])
        state = (
            engine.load(checkpoint, menu)
            if checkpoint.exists()
            else engine.initial(
                alphabet.SEEDS[alphabet.ORDER[pass_number]]["central"], menu
            )
        )
        journal = directory / f"candidates_{pass_number}.jsonl"
        # Equality is the complete canonical code; the hash is only an identifier.
        known = {
            growth.canonical_code(r["faces"]) for r in rows(journal, repair_tail=True)
        }

        def sink(faces, code):
            if code not in known:
                append(
                    journal,
                    dict(
                        candidate_id=hashlib.sha256(repr(code).encode()).hexdigest(),
                        faces=code[1],
                    ),
                )
                known.add(code)

        if not state["exhausted"]:
            engine.run(state, seconds, sink, checkpoint, stopped)
        result = summary(directory, pass_number, state, len(known))
        status = "ENUMERATION EXHAUSTED" if state["exhausted"] else "PAUSED"
        print(
            f"pass {pass_number}: {status}; candidates={len(known)}; "
            f"outer={state['stats']['outer_states']}; search={state['stats']['t']:.2f}s",
            flush=True,
        )
        return result["enumeration_exhausted"]


def parallel(directory, workers, seconds, stopped):
    """Independent processes own different passes; every child checkpoints on stop."""
    pending = iter(range(len(alphabet.SIGS)))
    jobs = []
    exhausted = False
    try:
        while (not exhausted or jobs) and not stopped():
            while not exhausted and len(jobs) < workers:
                number = next(pending, None)
                if number is None:
                    exhausted = True
                    break
                command = [
                    sys.executable,
                    str(HERE / "census.py"),
                    "exhaust",
                    str(number),
                    str(seconds),
                    "--cache-dir",
                    str(directory),
                ]
                jobs.append(subprocess.Popen(command))
            for job in jobs[:]:
                status = job.poll()
                if status is not None:
                    jobs.remove(job)
                    if status:
                        raise RuntimeError(
                            f"Census worker failed with exit code {status}"
                        )
            if jobs:
                time.sleep(0.2)
    finally:
        for job in jobs:
            if job.poll() is None:
                job.terminate()
        for job in jobs:
            job.wait()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=["grow", "exhaust", "all", "parallel", "status", "merge", "passes"],
    )
    parser.add_argument(
        "numbers",
        nargs="*",
        type=float,
        help="pass and seconds; or workers and seconds; or seconds for all",
    )
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args(argv)
    maximum = {"grow": 2, "exhaust": 2, "all": 1, "parallel": 2}.get(args.command, 0)
    if len(args.numbers) > maximum or any(not math.isfinite(n) for n in args.numbers):
        parser.error("Unexpected or nonfinite numeric argument")
    numbers = args.numbers
    if args.command in ("grow", "exhaust"):
        if (
            not numbers
            or not numbers[0].is_integer()
            or not 0 <= numbers[0] < len(alphabet.SIGS)
        ):
            parser.error("An integer pass number from 0 to 62 is required")
    if (
        args.command == "parallel"
        and numbers
        and (not numbers[0].is_integer() or numbers[0] < 1)
    ):
        parser.error("Worker count must be a positive integer")
    seconds = (
        (numbers[0] if numbers else 60)
        if args.command == "all"
        else (numbers[1] if len(numbers) > 1 else 60)
    )
    if seconds <= 0:
        parser.error("Chunk seconds must be positive")
    if args.command == "passes":
        for number, (label, word) in enumerate(zip(alphabet.ORDER, alphabet.SIGS)):
            print(
                f"{number:2d}  {label:15s}  defect={360-sum(alphabet.ANG[n] for n in word):3d}°"
            )
        return
    directory = cache(args.cache_dir)
    if args.command in ("status", "merge"):
        reports = sorted(
            (json.loads(p.read_text()) for p in directory.glob("pass_*_summary.json")),
            key=lambda r: r["pass_number"],
        )
        if args.command == "merge":
            atomic_json(
                directory / "merged.json",
                dict(passes=reports, geometric_certification_complete=False),
            )
        for report in reports:
            print(
                report["pass_number"],
                report["link"],
                "exhausted" if report["enumeration_exhausted"] else "paused",
                "candidates",
                report["candidates"],
            )
        return
    stop = False

    def request_stop(*_):
        nonlocal stop
        stop = True

    previous = {
        sig: signal.signal(sig, request_stop) for sig in (signal.SIGINT, signal.SIGTERM)
    }
    try:
        if args.command == "parallel":
            parallel(
                directory, int(numbers[0]) if numbers else 4, seconds, lambda: stop
            )
        else:
            passes = (
                range(len(alphabet.SIGS))
                if args.command == "all"
                else [int(numbers[0])]
            )
            for number in passes:
                while not stop:
                    if (
                        grow(directory, number, seconds, lambda: stop)
                        or args.command == "grow"
                    ):
                        break
                if stop:
                    break
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


if __name__ == "__main__":
    main()
