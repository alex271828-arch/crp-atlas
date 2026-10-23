# BI certifier 4.0 — final validation

## Release results

- **54 regression and adversarial tests pass** with the rebuilt helper: `validation/final-tests.log`.
- **206/206 saved certificates independently replayed by the final checker**: 198 OPEN_FOR_ALGEBRAIC and 8 EXCLUDED_NONFLAT; all match the earlier results. Evidence: `validation/final-replay/summary.json` and individual reports.
- **206/206 fresh Census 1.16 inputs completed** with new numerical estimates, new certificate proposals and independent standard-library replay: 198 OPEN_FOR_ALGEBRAIC and the same 8 EXCLUDED_NONFLAT. No unresolved cases. Evidence: `validation/final-census/completion.json`, `summary.json`, per-case certificates and `verified.json`.
- Exact canonical face-complex comparison matches every fresh input and verdict to the reference corpus, with no duplicate input complexes: `validation/reconciliation.json`.
- Full symmetry permutations and vertex orbits match the independent v3 graph matcher on all 206 inputs: `validation/symmetry-comparison.json`.
- The 1918 tracked original v3 files remain unchanged. The final validation manifest confirms that Python sources and the numerical binary were unchanged during the completed final run.

Saved-certificate replay was completed across two final-code passes: the first 15 reports, then the remaining 191 with a 900-second allowance. The earlier short-budget replay was stopped after large cases exceeded its five-minute harness limit under parallel load. Final reports are consolidated in `final-replay`; the remaining-case pass also retains its own summary in `replay-completion`.

The fresh validation uses four workers, 10 displayed digits, the unchanged gradient settings, and **600 seconds per stage**. The normal CLI default remains 300 seconds per stage. The longest complete candidate took 631.995 seconds across its stages. This is not a guarantee for other machines or inputs. Enumeration and algebraic realization were not rerun.

## Audit and tests

The tests cover exact outward operations against rational endpoint extrema at several precisions; the exact integer fast path; dot-product enclosure equivalence; cache behavior when precision changes; malformed rational inputs, labels and triangle lists; corrupted centers, inverse matrices and ridges; nonflat witnesses; improper coplanar inputs; full Catalan enumeration controls; hull metric transitions; preservation of unresolved statuses; timeouts; and proof checks under optimized Python. Numerical tests independently compare analytic gradients with finite differences and exercise symmetry, restarts and failure flags.

Both proof routes remain necessary. Full convex metric certification never asserts exact original-face flatness or CRP° membership. See `MATHEMATICAL_AUDIT.md` for the reviewed argument and its limitations.

## Measured improvements and code size

Symmetry enumeration over the 206-case comparison took 5.916 seconds with v3 and 0.418 seconds with direct face propagation, about 14.2× faster. An isolated interval dot-product benchmark measured 2.75× speedup for integer accumulation over separately rounded products. These are component measurements, **not an overall census speedup claim**. The integer endpoint optimization was motivated by the recorded large-case profile; it avoids repeated general rational construction without changing the arithmetic result.

Top-level Python executable statements: 1347 → 1289. Physical lines: 1509 → 2099; the increase comes from readable formatting and explicit checks. C++ gradient/helper sources compare identically after formatting normalization; they were rebuilt locally. NetworkX is removed from runtime dependencies.

## Development records

`final-replay`, `final-census`, `final-tests.log` and `reconciliation.json` are the release evidence. Other validation folders/logs are preliminary development evidence. A pre-optimization checker already replayed all 206 saved certificates, but a preliminary fresh batch hit five-minute proposal timeouts on large decagonal patches. Those cases were retained as unresolved, never excluded. Profiling identified repeated integer-to-Fraction conversion as a substantial cost; the exact integer path was then added, tested and subjected to the final complete replay and fresh run. Interrupted attempts remain preserved. Earlier compiler setup failures are retained in preliminary logs; `build-final.log` records the successful rebuild.

The release ZIP contains code, dependencies, examples, tests and validation summaries. Full local per-case proof artifacts remain in this release directory's `validation/final-census` and `validation/final-replay` directories.
