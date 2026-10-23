# Census v1.16 validation — release snapshot

Release date: 25 September 2026. At the user's request, the release is completed
without waiting for the slowest full-census passes. Those searches continue
in the background on the Mac. **The full 63-pass comparison is still pending.**

## Completed checks

- **45 regression tests passed**, in 42.434 seconds. This includes the retained
  growth/pruning regressions, all 120 catalogue controls, independent cyclic
  support and polygon-run enumerators, exact integer arithmetic, pinched
  fragments, canonicalization, journal recovery, interrupted expansion,
  checkpoint resume, locking, cache equivalence and CLI input validation.
- Four complete restricted-menu enumerations, independent of catalogue data,
  recover exactly the octagonal/decagonal prisms and antiprisms.
- Census executes with Python's site/package loading disabled, with no NumPy,
  SciPy, NetworkX, mpmath or historical realizer imported.
- All **717** original v1.15 source/result hashes checked during development
  remain unchanged.
- The real parallel run was interrupted and resumed from its own checkpoints.

## Full-census snapshot

**57/63 passes exhausted and matched v1.15 exactly**, with no completed-pass
mismatches. The currently saved journals contain 191 candidates;
none is outside the old candidate set. Unfinished passes: 41, 48, 49, 60, 61, 62.
The complete v1.15 reference contains 206 candidates.

These are snapshot counts, not a claim that the new census has completed or
that a pending pass will necessarily match. Live checkpoints and journals are
in `validation/full-census/`; the live process log is `validation/full-census.log`.
Read the per-pass done files for current completion, not this frozen snapshot.

When the remaining passes finish, update the comparison with:

```sh
python3.10 tools/compare_census.py \
  /Users/ferleger/Downloads/crp3456810_census_v1_15/cache_links_exact_extended \
  validation/full-census --output validation/census-comparison.json
```

The command exits successfully only if all 63 passes and all candidate sets
match, with no duplicate new records. It never reruns enumeration or BI.
Matching v1.15 is regression evidence, not a proof that both versions are free
of a shared mathematical error.

## Measured improvements

- Top-level runtime code: **3,335 → 1,855 lines** (44% fewer), despite expanding
  compressed statements into ordinary readable formatting. Runtime statements:
  2,980 → 1,246. This comparison includes the removed legacy/numerical modules.
- Canonicalization on 1,030 identical saved inputs: 2.957 → 2.198 seconds,
  approximately **26% less time**, with identical full codes.
- A fixed 12-fragment, box-limited geometry workload: approximately 3.21 → 2.48
  seconds, approximately **23% less time**, with identical interval covers.
  This isolates exact-endpoint sine caching; it is not an end-to-end census
  speedup claim. Python 3.14 also passed the same cover-equivalence benchmark.

## Evidence and limitations

`tests.log`, `audit-tests.log`, `restricted-census.json`, `standalone.json`,
`code-size.json`, `canonical-v15.json`, `canonical-v16.json`,
`sine-benchmark-final.json`, `environment.json`, and `original-preservation.json`
are under `validation/`. The benchmark fixture and tools allow local replay.
`pre-cache-census/` is an interrupted preliminary run, not the final-code census.

No full BI or algebraic rerun was needed or started. Necessary interval pruning
is retained, but candidate existence and original-face flatness remain separate
certification obligations. See `MATHEMATICAL_AUDIT.md` for the reviewed coverage
argument and its implementation-level, non-formal proof scope.
