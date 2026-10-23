# CRP° metric converter — fourth component

Convert saved exact algebraic realization data into explicit vertex coordinates in both a general positive-definite metric and an ordinary orthonormal Cartesian frame. Every successful output requires a separate coordinate-checking process. Existing Census, BI and algebraic packages are read-only inputs.

## Pipeline

1. Census 1.16 enumerates face complexes.
2. BI certifier 4.0 certifies convex metric realizations and excludes nonflat original faces.
3. Algebraic radii 3.0.11 and its saved earlier results supply exact radii and realization data.
4. This converter constructs explicit coordinates, then independently verifies the entire coordinate realization.

The input is a saved case directory containing `input.json` and `radii.json`. Radii alone do not specify the task: the converter also needs the labeled face incidence, assignment of vertices to spheres, primitive polynomial and its isolated real root. Source solver success flags and numerical coordinates are not proof premises.

## Use

```sh
python3.10 convert.py /absolute/path/to/saved-case --output /absolute/path/to/new-output
python3.10 check.py /absolute/path/to/new-output/coordinates.json
```

Install the tested dependency with `python3.10 -m pip install -r requirements.txt`. The output directory must be new. `--timeout` defaults to 600 seconds **per stage** (proposal and verification). A failed or interrupted stage is UNVERIFIED, never a nonexistence claim. The standalone checker can also be run under Python `-O`; proof conditions remain enabled.

## Coordinate convention

Vertex labels are those of the selected algebraic input. The catalog matches bodies to BI cases by canonical incidence; it does not relabel vertices into the BI input numbering.

`metric_coordinates[i]` is a three-vector qᵢ with squared distance `(qᵢ−qⱼ)ᵀG(qᵢ−qⱼ)`. `metric` is the exact symmetric positive-definite matrix G. The sphere center is recorded separately.

`cartesian_transform` is B with **BᵀB = G**. `cartesian_coordinates[i] = Bqᵢ` uses the ordinary Euclidean dot product. These are explicit algebraic Cartesian coordinates, not decimal approximations.

A sparse coordinate element maps a decimal bit mask to a base-field coefficient. Bit i contributes the **positive** square root rᵢ of `positive_radicals[i]`. For example, `{ "0": "1/2", "3": "alpha" }` means `1/2 + alpha*r0*r1`. The real value of alpha is fixed by the primitive polynomial and rational isolating interval. Coefficients may also contain explicit positive square roots of regular-polygon constants; the coefficient field is therefore not asserted to equal Q(alpha). Human-readable formulas are included and checked against this exact representation.

## Reconstruction

For radial input, an oriented fan triangulation of the original faces is developed from a first radial triangle. Three vertex rays form a nonorthogonal coordinate basis. Each new ray is obtained from exact scalar products with two already placed rays and a positive square root determining its oriented radial cone. All regular-polygon chords are derived symbolically. Missing consistency, nonpositive cones or unsupported algebra cause failure, not acceptance.

For an existing exact metric-coordinate certificate, its coordinates and metric are used only as a proposal. Both routes undergo the same independent global checks. Cartesian conversion uses a positive LDLᵀ factorization of the metric.

## Simplification and limits

The implementation reduces base coefficients exactly, recognizes equivalent radical square classes, shares positive radical definitions across all vertices and picks shorter factored coefficient expressions only after an exact equality check. This is a practical, exact simplification procedure. It does **not** claim the globally shortest expressions, a minimal realization field or a minimal Cartesian field degree. The supplied radius field can require extensions for coordinates; ordinary Cartesian coordinates can require radicals absent from a metric-coordinate description.

The source radial arithmetic module is vendored from algebraic radii 3.0.11. See `MATHEMATICS.md` for the verification argument, and `VALIDATION.md` for the completed release evidence. Python 3.10+ and SymPy are required. Acceptance uses exact equalities and rational interval signs. Numerical guesses inside a dependency’s field-embedding routine are only proposals: every embedded positive root must pass an exact square identity and positive-sign check.
