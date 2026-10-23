"""Isolate the interval dot-product change; this is not a census speed claim."""

import json
from pathlib import Path
import random
import statistics
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import exact_intervals as iv

rng = random.Random(460)
pairs = []
for _ in range(100):
    a = rng.randrange(-iv.SCALE, iv.SCALE)
    b = rng.randrange(-iv.SCALE, iv.SCALE)
    pairs.append((iv.I(a, a + 100, raw=True), iv.I(b, b + 100, raw=True)))


def old():
    return sum((a * b for a, b in pairs), iv.I(0))


def new():
    return iv.sum_products(pairs)


before, after = old(), new()
assert before.lo <= after.lo <= after.hi <= before.hi
samples = {"separately_rounded": [], "integer_accumulation": []}
for trial in range(4):
    methods = [("separately_rounded", old), ("integer_accumulation", new)]
    if trial % 2:
        methods.reverse()
    for label, function in methods:
        start = time.perf_counter()
        for _ in range(1000):
            function()
        samples[label].append(time.perf_counter() - start)
medians = {key: statistics.median(values) for key, values in samples.items()}
report = dict(
    dot_length=100,
    repetitions_per_trial=1000,
    bits=iv.BITS,
    samples=samples,
    median_seconds=medians,
    speedup=medians["separately_rounded"] / medians["integer_accumulation"],
    enclosure_not_wider=True,
    scope="Isolated interval dot products, not overall certification or census runtime.",
)
(ROOT / "validation/dot-benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
