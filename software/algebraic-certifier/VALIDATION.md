# Validation of algebraic radii solver v3.0.11

## Bivariate radial closure: body 129 (2026-09-17)

All **244 tests passed** in 41.191 seconds; see `TEST_RESULTS_v3.0.11.txt`.
Six new tests cover exact bivariate closure identities, unpaired denominators,
multiplicity bounds, the independently derived body-129 radial equations,
affine reduction from three radii to two parameters, and signed angle
identities with positive variable denominator factors. Existing deferred-pass
tests retain coverage of repeated reduction and timeout checkpoints.

Three-angle relations of total absolute multiplicity at most six are now
admitted with up to two radius parameters when the radial edge factors give a
provable square root of the denominator product. The existing cosine Gram
identity then avoids an additional squaring step. Unsupported relations are
retried after earlier exact radius reductions leave at most two parameters.
Larger or unpaired bivariate expressions remain deferred. No body identifier,
stored answer, coefficient fitting, geometry-extractor change, or extra time
allowance is involved. The compact-radical ceiling remains 12 seconds.

The original saved input for body 129 has two radius classes and a remaining
vertex relation of the form `3*A + B + C = 2*pi`. The old one-parameter gate
skipped its small polynomial consequence. The new consequence and the existing
radial relation isolate the two squared radii directly.

Fresh runs with the unchanged 120-second total limit:

| Body and mode | Solve seconds | Independent result |
|---|---:|---|
| body129-radius-only | 2.568 | CERTIFIED |
| body129-default | 5.703 | CERTIFIED |
| body153-radius-only | 7.540 | CERTIFIED |
| body163-radius-only | 27.860 | CERTIFIED |

Both body-129 runs use radial certification. The radius-only run enters no
geometry route; the default run uses the radial fallback after the geometry
extractor cannot complete its quadratic back-substitution. That separate
geometry limitation is unchanged.

Body 129 has two squared radii near 1.793347948 and 1.892133898. Their difference
d is near 0.098785950 and generates the degree-eight field with polynomial
`32*d^8 + 128*d^5 - 48*d^4 - 224*d^3 - 80*d^2 + 1`.
The readable report gives radii `a-d`, `a`, where
`a=(48*d^4+232*d^3+146*d^2-244*d+141)/(736*d^3+192*d^2-1376*d+196)`.
Exact presentation equivalence and nonzero denominators are independently
checked. The radial certificate verifies all 24 closed vertex fans with
winding number one, 34 positive face cones, 40 face-flatness checks, and
56 convex hinges.

The two previously difficult radial cases 153 and 163 still return SOLVED
and independently CERTIFIED, including their readable presentations.
Complete fresh outputs are in `VALIDATION/v3.0.11/`.

This release is installed in both `crp_algebraic_radii_v3` and the versioned
`crp_algebraic_radii_v3.0.11` folder. Earlier outputs and replaced files are
preserved. The full 198-body census was not rerun for this change.

## Exact coefficient arithmetic and body 111 (2026-09-17)

All **238 tests passed** in 39.877 seconds; see `TEST_RESULTS_v3.0.10.txt`.
Fourteen helper and package-import tests cover exact inverses, canonical
coefficients, algebraic zeros, rational and multiquadratic fields, explicit
polynomial variables, invalid inputs, and both import modes. Four budget tests
cover interior and noninterior centers, elapsed setup time, and the unchanged
overall wall-clock limit.

Constant algebraic coordinate matrices now use an explicit exact coefficient
field for inversion. Coefficients are reduced after changing Gram parameters,
so nested equivalent radical quotients do not accumulate in later rank
equations. The field is inferred from the exact coefficients; numerical
recognition and stored answer formulas are not used.

The geometry method receives the remaining search budget when the radial
fallback cannot apply. For interior centers, its old 65-second ceiling and
65% allocation remain. Time already spent is deducted and the existing
two-second reserve remains; the parent still enforces the requested total.

Body **111 returned SOLVED in 31.772 seconds** with
its original saved partition, default routes, and `--time-limit 120`. The
standalone checker returned **CERTIFIED**, including the readable presentation.
This already fits within the old 65-second geometry ceiling: arithmetic
normalization resolves the observed failure without relying on extra time.

Its field is Q(sqrt(5)), with squared radii
`1`, `(3+sqrt(5))/2`, and `(5+sqrt(5))/2`. The saved certificate uses
`4*alpha^2 - 2*alpha - 1`, alpha near -0.309017, and radii
`1`, `2-2*alpha`, `3-2*alpha`. The shorter presentation uses
`t^2 - 3*t + 1`, t near 2.618034, and radii `1`, `t`, `t+1`.

The independent exact geometry check confirms all 19 sphere memberships,
107 regular-face pair distances, 24 supporting face planes, and 21 maximal
facets, including exact tilings of the three subdivided coplanar facets.
It places the center on hexagonal face 15. The radial route's strictly
interior-center requirement and all certificate checks are unchanged.

Body **66 returned SOLVED in 6.184 seconds** using
the default corrected partition, followed by independent **CERTIFIED**.
All **10 portable examples returned SOLVED and independently CERTIFIED**,
with verified readable presentations. Their combined benchmark wall time,
including separate checks and presentation, was 196.298
seconds. Complete outputs are in `VALIDATION/v3.0.10/`.

The usual `crp_algebraic_radii_v3` installation and the versioned
`crp_algebraic_radii_v3.0.10` installation contain this release. Earlier
outputs and a backup of replaced files are preserved. The full 198-body
census was not rerun, and body 129 is not addressed by this validation.

## Default positive-radius partition for body 66 (2026-09-17)

All **220 tests passed** in 40.083 seconds; see `TEST_RESULTS_v3.0.9.txt`.
Thirteen new correction-selection tests cover the exact source match, census
aliases, unchanged custom partitions, stale geometry, modified replacement
metadata, and missing or empty catalogs.

The census launcher now applies a packaged correction to the original body-66
partition. The correction is matched by the original partition SHA-256 and
validated against the BI coordinates, face input, and replacement hashes. It
works both with the default collection and an explicit `--partitions` pointing
to that original collection. Exact solver and positivity checks are unchanged.

An audit of all 198 launcher selections found exactly one changed selection:
body 66, whose corrected groups have sizes 9, 4, and 2. Custom or modified
partition files remain as supplied. The original census files are preserved.

A fresh ordinary `run_body.py 66 --time-limit 120` run, with no custom partition
or route option, returned **SOLVED in 6.828 seconds**. The
independent checker returned **CERTIFIED**, including presentation equivalence.
The exact squared radii are alpha, alpha+1, and 5*alpha-1, where alpha is the
positive root near 0.904508 of `16*alpha^2 - 20*alpha + 5`.

Complete results, recheck, and the manual derivation are in
`VALIDATION/v3.0.9-body66/`. All v3.0.8 simplification improvements are retained.
The usual `crp_algebraic_radii_v3` installation and the versioned
`crp_algebraic_radii_v3.0.9` installation contain this release. The previous
usual installation was backed up; its historical output directories remain.
The full census was not rerun, and no changes for bodies 111 or 129 are claimed.

## Natural field generators and short rational expressions (2026-09-16)

All **207 tests passed** in 38.287 seconds; see `TEST_RESULTS_v3.0.8.txt`.
The new tests cover exact coefficient embeddings and real conjugates, proper
subfields and collective field generation, fallback on optional timeouts,
rational reconstruction, shared affine radius families, undefined denominators,
tampered presentations, and CLI budgets and forwarding.

The existing rational-coefficient root path is unchanged. A new optional
two-second attempt embeds the algebraic coefficient field into the isolated
closure-root field. It then proves, by closing a rational span under
multiplication, that the radii generate the entire field. Failure, a proper
subfield, or timeout retains the previous primitive-element construction.
No body identity, stored answer polynomial, or numerical recognition is used.

The bounded presentation stage retains short rational functions, tries simple
generators, and shares a common offset and affine coefficient relations.
Completed improvements are checked against the certified field, including
every denominator and the chosen real embedding. A presentation search timeout
does not change SOLVED. The independent certificate checker additionally checks
the presentation when present.

Fresh runs from the saved original inputs, with `--radius-only --time-limit 120`
and the default five-second presentation budget:

| Body | Solve seconds | Presentation seconds | Certificate display score | Readable score |
|---|---:|---:|---:|---:|
| 153 | 7.456 | 3.088 | 251 | 137 |
| 163 | 18.000 | 4.824 | 355 | 155 |

Both results returned **SOLVED**, and both standalone checks returned
**CERTIFIED**, including exact presentation equivalence. Body 163's presentation
search reached its optional budget and retained a verified shorter checkpoint.
The display score is the character count of the defining polynomial, shared
definitions, and radius formulas; conversion maps belong to the proof data.

Both certificate fields now use
`25*d^8 + 150*d^6 - 110*d^4 + 10*d^2 + 1`, with d=s1-s0.
The body-153 report gives a shared offset a and radii
`a`, `a+d`, `a+d*(-5*d^4-8*d^2+5)/(14*d^2-8)`.
Its old v3.0.7 certificate display score was 1239; the new readable score is 137.
The body-163 report finds an even shorter alternate generator t=s2-s1 with
polynomial `25*t^8 + 50*t^6 - 35*t^4 + 1` and radii
`a`, `a+b*t`, `a+(b+1)*t`, `a+(b+2)*t+1`.

The reports and complete saved results are in `VALIDATION/v3.0.8/body153/` and
`VALIDATION/v3.0.8/body163/`. This update is installed separately as
`crp_algebraic_radii_v3.0.8` because a v3.0.7 census run was active. The prior
package and its ongoing batch were not modified. The full census was not rerun.

## Retained partial reductions: body 163 (2026-09-16)

Completed deferred-angle reductions are now retained even when they do not
isolate a point. Unsupported fans are retried only after a strict decrease in
the number of radius variables. Previously materialized equations are carried
forward without regeneration; all passes share the same phase timer. A timeout
preserves the last completed reduction for the compact-radical fallback.

All **169 tests passed** in 33.544 seconds; see `TEST_RESULTS_v3.0.7.txt`.
Five new tests cover repeated materialization after dimension reduction,
retained progress in the fallback, termination without progress, and timeout
safety before and after a completed reduction. Three of these tests reproduce
the defects against v3.0.6 and pass with this repair.

Body 163 returned **SOLVED in 34.694 seconds** from its original saved input and
partition with `--radius-only --time-limit 120`. Deferred closure took 11.226
seconds. The trace records two passes: the first reduces two parameters to one,
and the second derives the missing univariate closure. No Gram solver is used,
and the compact-radical ceiling remains 12 seconds.

The standalone checker returned **CERTIFIED** for the saved body-163 result.
As a regression check, body 153 was rerun with the same options: **SOLVED in
18.190 seconds**, followed by standalone **CERTIFIED**.

Its squared radii are approximately 2.803966136816338, 3.537564845269763,
3.990953781197005, and 5.444342717124247, with respective multiplicities
20, 10, 10, and 2. Their common field is Q(d), of degree eight, where d is the
root in (733/1000,734/1000) of
`25*d^8 + 150*d^6 - 110*d^4 + 10*d^2 + 1`.

An independent manual derivation from the radius-angle equations gives the
same values and passed direct and standalone verification: 42 closed vertex
fans, 60 positive face cones, and 100 convex hinges. The derivation and fresh
solver outputs are in `VALIDATION/v3.0.7-body163/`. Existing batch results were
preserved; the full census was not rerun.

## Deferred one-parameter closure: body 153 (2026-09-16)

The solver now substitutes the derived one-parameter radius family into a
deferred three-angle relation before applying its bounded cosine Gram
construction. This handles total angle multiplicity up to six in one parameter
while retaining the existing multivariable expansion limits and independent
exact verification.

All **164 tests passed** in 33.972 seconds; see `TEST_RESULTS_v3.0.6.txt`.
Five new tests cover signed constant and univariate angle identities,
multivariable construction limits, the exact body-153 closure factor with no
denominator pole, and automatic substitution of a derived radius family in a
synthetic three-radius problem.

A fresh body-153 run using its original saved partition, `--radius-only`, and
the unchanged 120-second total limit returned **SOLVED in 18.516 seconds**.
Its missing closure took 1.083 seconds. The standalone saved-certificate checker
also returned **CERTIFIED**. No Gram search was performed and the 12-second
compact-radical limit was unchanged.

The three squared radii are approximately 2.803966136816338,
3.537564845269763, and 3.990953781197005, with group sizes 20, 10, and 10.
Their common field is Q(d), where d is the root in (733/1000,734/1000) of the
irreducible polynomial `25*d^8 + 150*d^6 - 110*d^4 + 10*d^2 + 1`.
This polynomial was derived from the original radius-angle equations.

An independent regular-polygon layer construction produces the same values and
checks all 240 within-face vertex-pair distances and all 40 sphere memberships
exactly. The full radius verifier additionally checks positivity, all 40 fan
closures, and all 90 convex hinges. Exact formulas, both derivations, and saved
certificates are included in `VALIDATION/v3.0.6-body153/`.

Existing batch results were preserved; the full census was not rerun.

## Radius-only mode and closure budget (2026-09-16)

Version 3.0.5 adds `--radius-only` to the solver and both launchers. It skips
Gram construction and solving while retaining exact candidate verification.
The existing strictly-interior center condition is enforced explicitly in this
mode. The compact signed-radical closure phase now has a 12-second ceiling,
subject to the remaining reduction budget and hard total wall-clock limit.

All **159 tests passed** in 33.603 seconds; see `TEST_RESULTS_v3.0.5.txt`.
Nine new tests cover route selection, rejected candidates, the interior-center
guard, equation generation, absence of Gram imports, and flag forwarding through
the solver parent process, census launcher, and benchmark launcher.

Fresh runs with `--radius-only --time-limit 120` produced:

| Census body | Result | Seconds | Route and independent check |
|---|---|---:|---|
| 68 | SOLVED | 7.657 | `radial_only`; standalone `CERTIFIED` |
| 153 | UNRESOLVED | 13.601 | Compact closure exhausted its 12-second allowance |

Neither run entered Gram construction or wrote a Gram derivation. Body 153's
trace records `budget_seconds: 12.0`, and no candidate reached verification.
The outputs are in `VALIDATION/v3.0.5-radius-only/`. Existing batch results were
preserved.

## Algebraic constant repair (2026-09-16)

The reducer now cancels constant residual equations and checks remaining
algebraic zero identities exactly. Both terminal solution paths use these
checks. Radical reduction treats constant numerators and denominators as
coefficients and factors nonconstant expressions with explicit generators.

Twelve focused regression tests cover algebraic denominator cancellation,
nested radicals and isolated-root identities, nonzero residuals as small as
10^-300, solved and contradictory linear systems, constant radicands, and
symbolic radicands with algebraic coefficients and denominators. They use no
numerical zero threshold.

All **150 tests passed** in 33.909 seconds; the complete output is in
`TEST_RESULTS_v3.0.4.txt`.

The seven bodies that previously reported constant-factorization errors were
rerun from their original saved inputs and partitions, with a 120-second limit
per body. Every successful saved certificate was checked independently.

| Census body | Result | Seconds | Independent recheck |
|---|---|---:|---|
| 67 | SOLVED | 6.324 | CERTIFIED |
| 68 | SOLVED | 72.824 | CERTIFIED |
| 117 | SOLVED | 6.281 | CERTIFIED |
| 153 | UNRESOLVED | 7.593 | — |
| 154 | UNRESOLVED | 34.830 | — |
| 163 | UNRESOLVED | 12.354 | — |
| 175 | UNRESOLVED | 74.790 | — |

The four unresolved cases now report that no isolated witnessed radius point
was extracted within supported elimination limits. None of the seven runs
reported the previous factorization exception or a traceback. The three solved
cases have squared radii 3/4, 7/4, and (5+2*sqrt(5))/4.

These fresh outputs and their rechecks are distributed in
`VALIDATION/v3.0.4-constant-regressions/`. Existing batch outputs retain their
original version and verdicts; the full census was not rerun for this repair.

## Coplanar subdivision repair (2026-09-16)

The geometric verifier now proves exact tilings of maximal supporting facets
by the listed regular polygons. It retains exact convexity, face-distance,
topology, and sphere-membership checks. Census 37 is the regression case: two
coplanar triangle pairs tile two rhombi, and its sphere center is on a square
face. The squared radii are 1/2, 3/2, and 11/6.

Focused tests cover this case with independently reversed and reordered faces,
standalone saved-certificate verification, and corrupted geometry or radii.
Additional facet-tiling tests reject overlaps, missing boundary vertices,
multiple boundary components, and nonconvex boundaries.

All **138 tests passed** in 32.46 seconds; the full output is in
`TEST_RESULTS_v3.0.3.txt`. A fresh census-37 run returned `SOLVED` in 0.69 seconds,
and its standalone saved certificate returned `CERTIFIED`. The older J92
geometric certificate also passed the updated standalone checker.

After installation in the local package, the affected completed batch entries
37, 39, 43, 44, and 60 were rerun. All five returned `SOLVED`, and all five
saved certificates independently returned `CERTIFIED`. Their previous outputs
were backed up before replacement. The existing batch continues normally.

## Previous v3.0.2 benchmark

Measured on the development Mac on 2026-09-15, using Python 3.10.2,
SymPy 1.13.1, NumPy 2.2.4, and python-flint 0.9.0. Existing BI coordinates were
used throughout; neither BI nor the census was rerun.

All ten examples were solved and their saved certificates passed the standalone
checker. Times include the solver's exact verification and output, but exclude
the second, independent check. Each run had a hard 120-second limit.

| Body | Seconds | Primitive field degree | Verification method | Independent recheck |
|---|---:|---:|---|---|
| CRP87 | 1.03 | 1 | radial | passed |
| CRP140 | 2.34 | 2 | radial | passed |
| J20 | 3.70 | 4 | radial | passed |
| J84 | 1.55 | 3 | radial | passed |
| J85 | 10.11 | 6 | radial | passed |
| J86 | 5.53 | 4 | radial | passed |
| J87 | 24.47 | 8 | radial | passed |
| J88 | 43.34 | 16 | radial | passed |
| J89 | 28.73 | 10 | radial | passed |
| J92 | 9.90 | 2 | exact_gram_geometry | passed |

The complete benchmark, including separate rechecks, took 174.71
seconds. Full output is in `VALIDATION/`. All 120 automated tests pass; their output
is in `TEST_RESULTS.txt`.

## J20 and general number-field root extraction

The former root routine failed when `Poly.lift()` encountered mixed rational
and algebraic coefficients. The replacement uses the field norm to obtain a
rational polynomial, isolates its real roots, and verifies candidate roots by
exact substitution into the original polynomial. This last step rejects roots
from other coefficient embeddings. Rational-coefficient inputs retain their
direct path. No dependency upgrade or body-specific algebra is required.

J20 now solves through the local-radius route. Its squared radii satisfy

    s0^2 - 5*s0/2 + sqrt(5)/10 + 17/16 = 0
    s1 = 7/4 + sqrt(5)/2

The certified branch is s0 = 5/4 + sqrt((5-sqrt(5))/10). The primitive
polynomial is

    1280*alpha^4 - 6400*alpha^3 + 10720*alpha^2 - 6800*alpha + 1381.

Nine focused regression tests cover the J20 quadratic, rational inputs, mixed
coefficients, different quadratic embeddings, cubic and multiquadratic fields,
rational roots, repeated norms, extremely close wrong conjugates, and recoverable
failures of unsupported exact embedding checks. Numeric proximity only selects
candidates; it never replaces the exact embedding test or final geometric check.

## J92 repair

The saved two-sphere center is on the hexagonal face. Earlier v3.0.0 rejected
this before entering the Gram solver. V3.0.1 retains that center/partition,
skips unsigned radial equations, and derives exact metric coordinates.
A generic final affine reduction also avoids unnecessary polynomial lifting
over the face's quadratic coefficient field.

The independent certificate proves all 102 regular-face pair distances,
20 planar supporting facets, 18 distinct vertices and sphere memberships, and
exactly classifies the common center as boundary. Only the hexagonal face (index 11 in input.json) contains
the center; the other 19 face-plane tests are strictly interior. The squared radii
are 1 and (3+sqrt(5))/2, so

Q(T) = T^2 - (5+sqrt(5))/2 T + (3+sqrt(5))/2.

The standalone checker rechecks both the exact geometry and Q without numerical
coordinates or the solver's derivation. Tests also reject incorrect radii,
modified coordinates, invalid metrics, nonconvex regular-face data, unsafe
expressions and incorrect Q. Boundary and exterior centers are tested.

## Scope

The geometric checker accepts coplanar subdivisions only after certifying a
complete convex facet tiling. Existing interior radial certificates continue
to support coplanar facet subdivisions as well.
The full 198-body collection was not rerun. Time limits or unsupported cases
produce UNRESOLVED and never establish non-realizability.

## Reproduce

```sh
/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 benchmark.py \
  --out out/worked-examples --time-limit 120 --overwrite --recheck
/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 \
  -m unittest discover -p 'test_*.py'
```

Use the supplied requirements with a suitable Python environment elsewhere.
