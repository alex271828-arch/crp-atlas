#!/usr/bin/env python3
"""Complete small homogeneous menus and compare to explicit prism constructions."""

import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from census_alphabet import fan
from link_engine import initial, run
from link_growth import canonical_code


def prism(n):
    return [tuple(reversed(range(n))), tuple(range(n, 2 * n))] + [
        (i, (i + 1) % n, (i + 1) % n + n, i + n) for i in range(n)
    ]


def antiprism(n):
    return [tuple(reversed(range(n))), tuple(range(n, 2 * n))] + [
        face
        for i in range(n)
        for face in ((i, (i + 1) % n, i + n), ((i + 1) % n, (i + 1) % n + n, i + n))
    ]


if __name__ == "__main__":
    result = []
    for n in (8, 10):
        for kind, constructor, word in [
            ("prism", prism, (4, 4, n)),
            ("antiprism", antiprism, (3, 3, 3, n)),
        ]:
            state = initial(fan(word)[0], (word,))
            found = set()
            started = time.monotonic()
            while not state["exhausted"]:
                run(state, 30, lambda faces, code: found.add(code))
            expected = canonical_code(constructor(n))
            record = dict(
                kind=kind,
                n=n,
                exhausted=True,
                candidates=len(found),
                match=found == {expected},
                seconds=time.monotonic() - started,
            )
            result.append(record)
            print(json.dumps(record), flush=True)
    (ROOT / "validation/restricted-census.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    raise SystemExit(0 if all(r["match"] for r in result) else 1)
