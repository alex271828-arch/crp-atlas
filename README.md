# CRP° Atlas and Certification

[**Explore the interactive Atlas**](https://alex271828-arch.github.io/crp-atlas/) · [Download this project](https://github.com/alex271828-arch/crp-atlas/archive/refs/heads/main.zip)

An interactive catalogue of 198 certified convex regular-faced polyhedra from the finite census, together with exact parametric prism and antiprism families. Every census body includes algebraic radii, explicit Cartesian and metric coordinates, face cycles, and a replayable exact coordinate certificate. The family views overlap some census entries and are not additional disjoint counts.

## Browse

| Atlas bodies | Category |
|---|---|
| 1–5 | Platonic |
| 6–18 | Archimedean |
| 19–23 | Antiprisms |
| 24–28 | Prisms |
| 29–120 | Johnson J1–J92 |
| 121–198 | Others |

Drag to rotate, scroll to zoom, and use the information tabs for exact expressions. Colors distinguish face sizes. Browser drawings use decimal approximations; certification comes from exact arithmetic. The online site loads each body's data when selected.

## Four software components

| Directory | Component |
|---|---|
| `software/census` | Census 1.16: fragment growth and mathematically sound filtering |
| `software/bi-certifier` | BI 4.0: interval certificates and nonflatness exclusion |
| `software/algebraic-certifier` | Algebraic radii 3.0.11: exact realization machinery |
| `software/metric-converter` | Metric converter 1.0: explicit metric and Cartesian coordinates and their exact checker |

The mathematical kernels are exported from the existing releases. Build products, private working directories, interrupted runs, caches and platform binaries are omitted. Third-party source needed for the BI build is included. Historical package documents describe their own release-time snapshots; the final coordinate ledger and census comparison are the evidence for this publication.

The fresh publication replay passed **198/198 coordinate certificates and 8/8 exclusions**. See [the replay report](validation/publication-replay.json). All 201 hosted website files were downloaded and matched against the release.

## Verify the published certificates

Python 3.10 or later is required. From the repository root:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r software/metric-converter/requirements.txt
python3 verify.py --jobs 4
```

This independently checks all 198 coordinate realizations and the eight BI exclusions. It does not trust stored success flags. The coordinate checker verifies regular faces, support and facet tilings, sphere memberships, the Cartesian isometry and printed formulas. A timeout or any failed proof condition is reported as unverified, never as nonexistence. Expect several minutes or longer depending on hardware.

For a quick check of one Atlas body:

```sh
python3 verify.py --body 1
```

For direct verification:

```sh
python3 software/metric-converter/check.py certificates/001/coordinates.json
```

All published body numbers refer to the Atlas order above. Historical pipeline identifiers survive only inside provenance records. `catalogue.json` binds the public numbering to certificate hashes and exact reference-complex matches.

## Reproduce and inspect

- **Census:** standard-library Python. Run `python3 software/census/census.py passes`; see its README for the full exhaustive run. Completed candidate journals and the final comparison are under `data/census`.
- **BI:** C++17 and CMake are needed for the proposal engine; see `software/bi-certifier/README.md`. The independent certificate checker uses Python's standard library. Eight excluded candidates and their selected certificates are under `exclusions`.
- **Algebraic realization:** see `software/algebraic-certifier/README.md` and `solve.py --help`. A new solver run requires face incidence, BI numerical guides and a sphere partition. Numerical guides are proposals, not proof premises. This initial publication includes the saved exact inputs and radii used by the coordinate converter under `data/algebraic-inputs`, rather than every historical numerical search artifact.
- **Coordinate conversion:** `python3 software/metric-converter/convert.py data/algebraic-inputs/001 --output /path/to/new-output` reconstructs and verifies a new coordinate certificate.

The full software is provided; rerunning the entire census is separate from replaying these finite certificates. Historical absolute paths in source-provenance fields are descriptive and are not needed by the independent coordinate checker. The algebraic convenience selector requires an explicit `--partitions` directory in this public edition.

## Mathematical scope

The finite census uses face sizes 3,4,5,6,8,10, allows coplanar prescribed faces, requires positive vertex curvature and proper face intersections, and identifies reflections. The project proofs for sizes 7, 9 and ≥11 are under `scope-proofs`. The uniform algebraic family construction and supporting-plane argument are in [FAMILIES.md](FAMILIES.md). These documents and the census's exhaustiveness argument are distinct from the per-body coordinate certificates.

The expressions are exact with specified real embeddings and positive radicals. No globally shortest-expression or minimal-field claim is made.

## Local website

```sh
python3 -m http.server 8000 --directory docs
```

Open `http://localhost:8000`. The site needs HTTP because it loads the selected body's JSON. No server-side computation, account or installation is needed by public website visitors.

## Attribution and licensing

See [LICENSES.md](LICENSES.md). This publication preserves third-party licenses and does not assign a new blanket license to the research code. The reference catalogue has its own noncommercial-use terms. A repository-wide open-source license for original project code has not yet been selected.
