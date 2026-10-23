# BI certifier 4.0

Given a closed regular-polygon face complex, propose its convex intrinsic-metric realization numerically, then independently verify a rational interval certificate. Certify nonflat original-face defects when possible and retain every other candidate for algebraic realization.

This release is a separate cleanup of `crp_bi_numeric_v3`. Previous software and results are preserved. The numerical search still uses full face-complex vertex orbits and squared-radius gradient descent. The certifier uses no symmetry assumptions, algebraic coordinate recognition, stored body formulas, or catalogue acceptance decisions.

## Run on this Mac

```sh
cd /Users/ferleger/Downloads/crp_bi_certifier_v4
python3.10 run.py --input examples/cube.json --all
python3.10 run.py --census /Users/ferleger/Downloads/crp3456810_census_v1_16 --list
python3.10 run.py --census /Users/ferleger/Downloads/crp3456810_census_v1_16 --all
```

`--census` accepts a directory containing `candidates_*.jsonl`, a Census 1.16 installation with `validation/full-census`, or the older `cache_links_exact_extended` layout. Only candidate files supply input; accepted/rejected catalogues are never consulted. Use the identifier printed by `--list` with `--id`. Candidate display hashes and historical catalogue numbers are different identifier systems.

`--input` accepts a face-cycle list, an object containing `faces`, a list of such objects, or JSONL. Vertices must be contiguous integers starting at zero. Supported unit regular face sizes are 3,4,5,6,8,10. The closed sphere, positive-curvature and face-intersection conditions are verified. Coplanar neighboring faces are allowed.

Every run creates a new output directory. `--output-dir` chooses its name but cannot overwrite an existing directory. No run resumes by silently trusting an old result. `--timeout` is a positive finite number of seconds **per stage per candidate**, default 300. Its exhaustion never excludes a candidate.

`--digits` (3–500, default 10) controls outward decimal display and Cartesian refinement precision. Higher requested display precision does not guarantee arbitrarily narrow hull certificates: the hull producer proposes double-precision centers. Exact rational interval endpoints are always saved. Numerical controls and the explicitly uncertified `--numerical-only` diagnostic are listed by `--help`.

## Verdicts

| Status | Meaning |
|---|---|
| `EXCLUDED_NONFLAT` | A convex metric realization is certified and an original-face chord or height defect excludes zero. |
| `OPEN_FOR_ALGEBRAIC` | A convex metric realization is certified; exact flatness of the original faces remains unproved. |
| `UNRESOLVED_CONVEXITY` | A metric-equation root is enclosed but its convex identification remains unproved. |
| `UNRESOLVED_REALIZATION` | No completed endpoint proof; keep the candidate. |
| `NUMERICAL_ONLY` / `UNRESOLVED_NUMERICAL` | Certification was explicitly skipped; no proof claim. |

No BI verdict asserts CRP° membership or exact original-face flatness. A zero-containing interval is not an equality proof. Numerical failure, a resource cap or an unsupported metric triangulation is not nonexistence.

## Files to use

`summary.json` lists results; `next-stage.jsonl` contains every nonexcluded candidate. Each case records the actual input, untrusted numerical estimate, selected certificate, independent `verified.json`, exact coordinate enclosures and original-face measurements. Read the selected certificate and `verified.json`, not an earlier unsuccessful proposal report. Numerical and certified coordinate frames are explicitly distinguished. `manifest.json` records source and binary hashes and run settings.

For compatibility, verified reports retain empty `internal_folds` and explanatory `internal_fold_note` fields. Version 4 uses original-face height/chord witnesses; it no longer computes optional fold angles. The handoff schema remains version 2 because its fields and meaning are unchanged; verification reports identify `crp-bi-certifier/4`.

Replay a certificate independently:

```sh
python3.10 -S check.py PATH_TO_CERTIFICATE --output NEW_REPORT.json
```

The checker needs only Python's standard library. Proof conditions are explicit checks and remain enabled under `-O`. It accepts the earlier `interval-root` and `hull-metric-interval` certificate formats. Supplied success flags are ignored.

## Code and proof map

| File | Responsibility |
|---|---|
| `exact_intervals.py`, `field.py` | Outward integer interval arithmetic; exact regular-polygon constants. |
| `metric.py` | Input topology, intrinsic metric and flip replay. |
| `proof.py` | Rational parsing, proof failures and the common contraction check. |
| `check.py` | Strict convex support, original-face witnesses and verdicts. |
| `hull_metric.py` | Coplanar-patch neighborhood and complete piecewise metric derivative bounds. |
| `produce.py`, `hull_proposal.py` | Untrusted numerical certificate proposals. |
| `symmetry.py` | Direct face-adjacency automorphism enumeration for numerical optimization. |
| `gradient_embedder.h`, `metric_probe.cpp` | BI geometry and squared-radius gradient search. |
| `run.py` | Isolated stages, independent replay, reports and handoff. |

Read [MATHEMATICAL_AUDIT.md](MATHEMATICAL_AUDIT.md), [HULL_METRIC_PROOF.md](HULL_METRIC_PROOF.md), [GRADIENT_METHOD.md](GRADIENT_METHOD.md), [CHANGELOG.md](CHANGELOG.md) and [VALIDATION.md](VALIDATION.md). They distinguish proof obligations, numerical heuristics, measured coverage and remaining limitations.

## Dependencies and rebuilding

Python 3.10 or newer; NumPy, SciPy and mpmath for proposals. NetworkX is no longer a runtime dependency. The included macOS binary is platform-specific. The bundled Geometry Central and Eigen sources retain their licenses.

```sh
python3.10 -m pip install -r requirements.txt
./build.sh
python3.10 -m unittest discover -s tests -v
```

`build.sh` uses CMake (or `CMAKE=/path/to/cmake`) and a C++17 compiler. On this Mac it can use the installed Command Line Tools if the selected Xcode cannot run; it does not change system settings or accept licenses. CMake includes libc++ headers from recent SDKs when present. `BUILD_JOBS` defaults to 2.

Full saved-certificate replay and fresh-census validation tools are under `tools/`; they create new result directories and leave their inputs intact. Their corpus/reference paths are evidence inputs only, never part of the production solver.
