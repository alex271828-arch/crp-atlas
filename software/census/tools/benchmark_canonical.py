#!/usr/bin/env python3
"""Time a canonicalizer on identical saved complexes, in an isolated process."""

import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, sys.argv[1])
from link_growth import canonical_code

faces = []
for file in sorted(Path(sys.argv[2]).glob("candidates_*.jsonl")):
    faces.extend(json.loads(line)["faces"] for line in file.read_text().splitlines())
started = time.perf_counter()
codes = [canonical_code(f) for _ in range(5) for f in faces]
print(
    json.dumps(
        dict(
            calls=len(codes),
            seconds=time.perf_counter() - started,
            digest=hashlib.sha256(repr(codes).encode()).hexdigest(),
        )
    )
)
