#!/usr/bin/env python3
"""Compare cached/uncached sine arithmetic on fixed, box-limited proof workloads."""

import json
import pickle
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import certified_geometry as C
import link_growth as G

samples = pickle.loads((ROOT / "validation/geometry-benchmark-input.pkl").read_bytes())
samples = [(faces, G.analyze(faces, menu), region) for faces, menu, region in samples]
cached = C.sine


def raw(c):
    return C.sqrt(1 - C.square(C.cosine(c)))


trials = []
reference = None
for trial in range(4):
    use_cache = bool(trial % 2)
    for name in ["_quad_cached", "_pent_cached", "quad_domain", "pins", "_sine_bounds"]:
        getattr(C, name).cache_clear()
    C.sine = cached if use_cache else raw
    started = time.perf_counter()
    results = [C.solve(f, i, r, time_cap=10000, max_boxes=25) for f, i, r in samples]
    elapsed = time.perf_counter() - started
    if reference is None:
        reference = results
    if results != reference:
        raise AssertionError("Cache changed a proof cover")
    trials.append(dict(cached=use_cache, seconds=elapsed))
report = dict(samples=len(samples), identical_covers=True, trials=trials)
(ROOT / "validation/sine-benchmark-final.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps(report, indent=2))
