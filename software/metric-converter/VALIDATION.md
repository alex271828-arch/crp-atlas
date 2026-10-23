# Version 1.0.0 — release validation

All **198 of 198** survivors have explicit exact metric and Cartesian vertex coordinates, independently verified by the frozen release checker.

- All 198 labeled inputs and algebraic data match their source files by content and SHA-256.
- Canonical face-incidence matching establishes a bijection with the 198 BI 4.0 open survivors. Historical algebraic numbers and BI run numbers are different; consult the ledger.
- Input selection: 195 successful saved 3.0.7 cases plus the later repairs of historical cases 066, 111 and 129. There are 179 radial proposals and 19 proposals from saved exact Gram geometry.
- The final replay used a fresh checker process per case, six workers, and unchanged source hashes throughout. Elapsed wall time: 10.67 minutes.
- A fresh end-to-end conversion of pentagonal case 064 passes with the strengthened arithmetic; see `validation/strengthened-end-to-end`.
- **31 tests pass**, including deliberately corrupted metrics, coordinates, Cartesian maps, sphere assignments, printed formulas, malformed arithmetic, imaginary coefficients, invalid root isolation, wrong field-embedding proposals and topology. Proof rejection also works under Python optimization.
- Environment: Python 3.10.2, SymPy 1.13.1, macOS ARM64. See `validation/environment.json`.

## Evidence

`validation/final-198/completion.json` is the aggregate completion record. `summary.json` gives per-case checks and timings; `manifest.json` fixes the checker hashes and original input hashes. `incidence-ledger.json` supplies the exact BI cross-reference. Each numbered subdirectory contains the self-contained `coordinates.json`, the separate `verification.json`, checker log, and readable `coordinates.txt`.

The initial development trial is retained locally under `validation/all-198-trial` but excluded from the release archive. Code changed during that preliminary trial. An intermediate replay also predates the additional exact field-embedding checks and is excluded. Release evidence is the subsequent full replay with frozen code, not the preliminary completion record.

## Scope of the guarantee

Acceptance proves an exact full-dimensional convex realization with the prescribed regular flat faces, sphere memberships, and a verified isometry between metric and ordinary Cartesian coordinates. No floating-point tolerance is used for acceptance. Unsupported arithmetic, timeouts or unresolved signs cannot produce acceptance.

This is an exact-arithmetic computational certificate, not a proof-assistant formalization. Its trusted base includes Python, SymPy, the documented arithmetic and the standalone checker. Algebraic field minimality and globally shortest coordinate expressions are not claimed. The earlier three software packages are unchanged.
