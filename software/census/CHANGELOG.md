# v1.16 changes from v1.15

## Removed

- Numerical final fitting/classification, approximate strict/coplanar decisions,
  cap descriptions, hull naming and their five secondary ledgers.
- Catalogue access in the runtime and its belt-search dependency chain:
  `belt_seed.py`, `belts_clean.py`, `crf_realize.py`, `final_validation.py`.
- Superseded simple-boundary topology, curvature/run and forced-face routines.
  Their run-DP reference survives only as an independent test oracle.
- Historical source copies, unused NetworkX import and never-used planarity
  statistics. Prior releases retain all historical material.
- Redundant CLI flags and the older four-face filename alias.

## Simplified without changing mathematical rejection conditions

- One exact curvature function with caching, shared by local and boundary tests.
- Forced-face closure runs to a fixed point; monotone face addition replaces
  the arbitrary 92-step cap. Only intersections involving the new face change.
- One edge lookup per canonicalization orientation instead of per anchor.
- Bounded exact-endpoint sine caching; the interval formula is unchanged.
- Ordinary formatted control flow and documented geometric formulas replace
  compressed multi-statement lines. Proof-bearing branches are retained.
- Candidate-only CLI and journal handling are separate from the growth engine.

## Operational repairs

- Expansion is computed before committing agenda and seen-set changes. A sink
  or generator exception cannot remove an unexpanded item from memory.
- All step outcomes, including rejections and duplicate states, reach the
  periodic checkpoint check.
- A crash-torn journal tail can be repaired under the pass lock and replayed.
  Complete malformed records cause an explicit error.
- Checkpoints and JSON snapshots are flushed/fsynced before atomic replacement.
- A shared configuration lock prevents racing parallel initializations.
- CLI numeric validation rejects NaN, infinity, fractional worker counts and
  extra arguments. Signal handlers are restored after programmatic CLI calls.

## Preserved

All 63 passes and their order, the finite six-size face alphabet, positive
vertex curvature, coplanar neighboring faces, injective all-face growth,
pinched links, coupled gap rows, exact canonical equality, conservative
interval rejection, parent covers, and independent per-pass execution.

The old completed census was not reused as a search oracle. Fresh validation
compares candidate complexes only after their new passes have exhausted.
