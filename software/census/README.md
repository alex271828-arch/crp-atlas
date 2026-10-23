# Census v1.16

A standalone, resumable fragment-growth enumerator for CRP° with regular
unit-edge faces of sizes **3, 4, 5, 6, 8, 10**. Coplanar neighboring faces are
included. Vertices have strictly positive curvature. Distinct faces meet in
at most one vertex or one common edge. Reflections are identified.

This is stage 1: its output consists of **candidate face complexes**, not a
list of geometrically certified solids. Pass the journals to the separate BI
and algebraic stages. The infinite prism/antiprism families outside the finite
face alphabet require their separate mathematical treatment.

## Validation status

45 tests pass. The completed full-census passes match v1.15; the slow remaining
passes are still running at release time. See `VALIDATION.md` for the dated
snapshot and the command that checks the eventual complete result.

## Run

Only the Python standard library is needed (Python 3.10 or later on macOS/Linux):

```sh
cd /Users/ferleger/Downloads/crp3456810_census_v1_16
python3.10 census.py passes
python3.10 census.py parallel 4 60
```

The latter exhausts all 63 passes with four processes. Each process saves
resumable checkpoints. `60` is the search allowance for a chunk, not a total
run limit; the current indivisible expansion can overrun it. Ctrl-C requests
a checkpoint after that expansion; stopping the parallel parent also stops
its children gracefully.

Other commands:

```sh
python3.10 census.py grow 41 60        # one chunk of pass 41
python3.10 census.py exhaust 41 60     # finish pass 41
python3.10 census.py all 60            # finish all passes sequentially
python3.10 census.py status
python3.10 census.py merge
```

Use `--cache-dir /absolute/path` for a separate run. The default is `cache/`
inside this release. A matching run resumes automatically. Source hashes,
version and alphabet prevent accidental reuse of incompatible checkpoints.
Never load untrusted pickle checkpoints.

`census3456810.py` is a filename-compatible entry point to the new CLI.
The obsolete `--crp`, `--cpr`, `--extended`, `--defer-validation`, and `validate`
options are removed. The algorithm always enumerates the full CRP° candidate
class, and final certification is always separate. No accepted/rejected
numerical ledgers or catalogue names are produced.

## Output

For each pass K:

- `candidates_K.jsonl`: canonical, consecutively labelled face cycles and a
  SHA-256 identifier. Deduplication compares complete exact codes, not hashes.
- `pass_K_state.pkl`: frontier, rooted intermediate states, interval covers,
  exact seen sets and accumulated statistics.
- `pass_K_summary.json`: progress and candidate count.
- `pass_K_done.json`: written when that pass's agenda is empty.

Use the cache directory itself as the BI input, for example
`--census /Users/ferleger/Downloads/crp3456810_census_v1_16/cache`.

The journal schema remains compatible with the BI census reader:
`candidate_id` and `faces`. Vertex numbering is canonical in the new journals;
compare complexes by canonical code rather than by literal old labels.
`geometric_certification_complete` remains false even after enumeration ends.

## Code map

| File | Responsibility |
|---|---|
| `face_types.py`, `census_alphabet.py` | Exact angle units, curvature, 63 ordered passes and seeds |
| `link_growth.py` | Partial vertex links, all legal single-face additions, exact canonical codes |
| `boundary_filter.py` | Face-intersection test and cyclic boundary support |
| `surface_pruning.py` | Unique rotations, coupled missing words, curvature and forced-face closure |
| `certified_geometry.py` | Outward integer intervals, spherical-link transfers, interval covers and convex development |
| `rim_geometry.py`, `geometry_plan.py` | Incomplete-link chord bounds and free-variable planning |
| `link_engine.py` | Transactional outer/inner growth and checkpoints |
| `storage.py`, `census.py` | Journals, locks, configuration and parallel CLI |

The old belt engines, numerical realizer, catalogue naming, duplicated
simple-boundary filters, dead statistics and legacy source tree are not part
of this release. Catalogue data exists only under `testdata/` for regression
checks; runtime code never reads it.

See `EXHAUSTIVENESS.md` for the coverage argument, `MATHEMATICAL_AUDIT.md` for
the soundness review, and `VALIDATION.md` for observed results and limitations.
The original v1.15 directory and its completed census are preserved.

## Verification

The engine itself needs no third-party packages. Regression tests additionally
use NumPy and mpmath, already installed in this Mac's Python 3.10 environment:

```sh
python3.10 -m unittest discover -p 'test_*.py' -v
python3.10 tools/compare_census.py \
  /Users/ferleger/Downloads/crp3456810_census_v1_15/cache_links_exact_extended \
  validation/full-census --output validation/census-comparison.json
```

The comparison checks every pass and every candidate, including candidates
later excluded by BI. Matching v1.15 is regression evidence; it is not an
independent proof of the older or newer engine's mathematical completeness.
