# Algebraic radii solver v3.0.11

This version derives exact algebraic squared radii from a closed census face
complex and saved BI numerical coordinates. It does not rerun BI or the census.
It contains no Johnson-body lookup, stored answer polynomial, or numerical
algebraic recognition.

## Run

On this Mac:

```sh
cd /Users/ferleger/Downloads/crp_algebraic_radii_v3
BODY=J89
python3 run_body.py "$BODY" --time-limit 120 --overwrite
```

Change only `BODY`: it accepts a Johnson label, the original census number, or
the full CRP identifier. Census `87` and Johnson `J87` are different bodies.
The selector reads the matching saved BI coordinates and sphere partition.
Outputs go to `out/<census-number>-<CRP-identifier>/` within this package.
Use `--show-command` to inspect the selected paths without starting a run.
Use `--partitions /path/to/saved/partitions` to override the local collection.

### Two-parameter radial angle closure (v3.0.11)

The deferred radial solver now handles three-angle relations with total
absolute multiplicity at most six and at most two radius parameters when
the radial edge factors give a provable square root of the denominator
product. This uses the existing cosine Gram identity without an additional
squaring step. Earlier affine radius reductions can also unlock this case.
Relations outside these bounds remain deferred; the 12-second compact-radical
limit and full exact geometric verification are unchanged.

This resolves body 129 from its original partition through the radial route:

```sh
python3 run_body.py 129 --radius-only --time-limit 120 --overwrite
```

The fresh run returned SOLVED in 2.57 seconds and passed the independent
checker. The default routes also solve it, using the radial fallback.
Its two squared radii are approximately 1.793347948 and 1.892133898,
and their difference generates a degree-eight field with polynomial
`32*d^8 + 128*d^5 - 48*d^4 - 224*d^3 - 80*d^2 + 1`.

All earlier fixes are retained. This release is also available in the separate
`crp_algebraic_radii_v3.0.11` directory.

### Exact coefficient arithmetic and geometry budgets (v3.0.10)

Constant algebraic matrices are inverted in their exact coefficient number
field. After changing Gram parameters, polynomial coefficients are reduced
in a shared exact field, with the unknown parameters kept as polynomial
variables. This prevents equivalent radical expressions from growing through
successive coordinate transformations. These are general arithmetic changes;
the solver still derives and independently checks every candidate.

When the supplied center is on the boundary or outside the body, the radial
fallback is unavailable. The geometry search now receives the remaining search
budget in that case. With an interior center, it retains the existing cap of
65 seconds and at most 65% of the remaining search budget, preserving time for
the radial fallback. The requested overall time limit remains enforced.

Body 111 now solves from its existing input and partition in about 32 seconds
with the default routes and a 120-second limit. Its independently certified
squared radii are `1`, `(3+sqrt(5))/2`, and `(5+sqrt(5))/2`.

```sh
python3 run_body.py 111 --time-limit 120 --overwrite
```

Body 66's positive-radius partition and the field and presentation improvements
are retained.

### Corrected numerical sphere guides (v3.0.9)

The census launcher automatically applies packaged corrections to known saved
partitions. This release includes body 66's verified positive-radius cover,
with three vertex groups of sizes 9, 4, and 2. Its original saved partition
placed the common center at vertex 7 and therefore had an exact zero radius.

```sh
python3 run_body.py 66 --time-limit 120 --overwrite
```

This works with the default collection and with `--partitions` pointing to the
original collection. The correction catalog matches the original partition's
SHA-256 hash and verifies the BI coordinates and face-file hashes before
selecting the replacement. Other partitions, including changed custom inputs
for the same body, are left as supplied. No original census file is overwritten.
The selected correction is announced by the launcher, and the solver records
the replacement file and hash in its source provenance. Its packaged JSON also
records the original partition hash.

The correction supplies only numerical guide data and vertex groups. Exact
radii are still derived and certified anew; zero-radius spheres remain
rejected. To reproduce the original input deliberately, pass it directly to
`run.sh --spheres /path/to/original.json` with the corresponding case and BI
arguments. The low-level solver honors that explicit input.

Add `--radius-only` to skip the geometry/Gram solver and go directly from the
short radius prepass to the full radius equations:

```sh
python3 run_body.py 153 --radius-only --time-limit 120 --overwrite
```

The same flag is accepted by `solve.py`, `run.sh`, and `benchmark.py`. Radius
equations require a strictly interior sphere center: with this flag, a supplied
boundary or exterior center returns `UNRESOLVED` with an explicit reason. No
geometry route is attempted, including with `--equations-only`. Successful
radius candidates still undergo the full exact verification.

The ten included examples are portable, requiring no external project files:

```sh
./run.sh --case examples/J89 --spheres examples/J89/spheres.json \
  --out out/J89 --time-limit 120 --overwrite
```

The examples are J20, J84 through J89, J92, CRP87, and CRP140. They contain face incidences
and numerical witnesses, with provenance, but no exact answers.

For another input, `--case` accepts a BI directory containing `input.json`
(faces) and `numerical.json` (positions), or a face JSON with a separate `--bi`.
Faces have consecutive vertex labels starting at zero. Allowed face sizes are
3, 4, 5, 6, 8, and 10. Without `--spheres`, a bounded numerical search proposes
a sphere partition. There is no imposed four-sphere limit.

Python 3.10+, SymPy and NumPy are required. `run.sh` selects the installed
Framework Python on this Mac; elsewhere install `requirements.txt` and use
Python with those dependencies. `CRP_RADII_PYTHON` overrides the wrapper's choice.

## Method

A short prepass first tries inexpensive radius identities. By default, the main path then
works with the metric before choosing the common sphere center.

1. Build the exact centroid Gram matrix from unit edges, regular-face diagonals,
   all combinatorial automorphisms, the centroid identity, and the affine
   dependencies within each regular polygon. Solve these linear constraints.
2. Split this affine Gram family using exact projectors for commuting
   involutions. Choose useful exact Gram entries as parameters. For a rank-r
   block, project onto r independent anchors and impose its Schur complement
   as zero. These are the distance-intersection equations in inner-product
   form: linear projection followed by residual norm equations.
3. Reduce the remaining local equations, preserving small degree where possible.
   Rank-one blocks use outer products; exact reduction recovers invariant inner
   products without retaining unnecessary coordinate sign choices. Higher-rank
   blocks seek constant local metrics among projected vertex differences, then
   solve line/conic intersections through a known point. A bounded small
   polynomial ideal and rational root isolation finish the remaining closure.
4. Recover the common sphere center from the linear equal-distance equations
   for the proposed vertex groups. If these leave an affine family, choose its
   center nearest the vertex centroid. Calculate the exact squared radii.
5. Independently verify the original body. Interior centers use the radial
   certificate. For a boundary or exterior center, a separate exact verifier
   checks metric positivity, every face distance and plane, convex supporting
   planes, and sphere membership. Only a successful check permits `SOLVED`.

A rejected algebraic candidate does not prevent trying the next route. A
bounded fallback uses exact local metric patches and short radial-angle
identities. Dense angle closures are retained as templates and instantiated
only after simpler relations have reduced the variables. V3 does not invoke
v2's special two-square coordinate template.

With `--radius-only`, the solver skips the metric construction and tries this
radial closure immediately after the prepass. Version 3.0.5 raises the compact
signed-radical closure ceiling from 6 to 12 seconds in both modes. The total
wall-clock limit and allocation from the remaining reduction budget can shorten
that phase.

Numerical coordinates guide the sphere partition, symmetry representation
ranks, nonsingular pivots, polynomial factors, and real embeddings. They never
provide exact coefficients. These choices make the search incomplete: a wrong
or unsupported branch can return `UNRESOLVED`, but cannot certify impossibility.
Exact positive-semidefinite rank-three checks validate a completed Gram matrix;
the independent verifier then checks the original body from its exact radii
or its exact metric coordinates and sphere center.

The Gram search has explicit size and time bounds. Bodies with larger or less
useful symmetry decompositions may require the radial fallback or remain
unresolved. This is a solver prototype, not a complete decision procedure for
all 198 survivors. Passing BI alone is not treated as an exact certificate of
regular flat faces.

## Deferred one-parameter angle closure

Version 3.0.7 preserves each completed partial reduction of the deferred angle
equations, even when one or more radius parameters remain. If a reduction
removes a variable, previously unsupported fans are retried with the smaller
radius family. Materialized equations are retained, and every pass shares the
same phase time budget. A timeout retains the last completed reduction for the
subsequent compact-radical route.

This resolves body 163 through `--radius-only`: one fan first reduces two
radius parameters to one, and a second fan then supplies the final closure.
The common squared-radius field is the same degree-eight field as body 153.
The 12-second compact-radical ceiling is unchanged.

Version 3.0.6 substitutes the already-derived radius family into a deferred
angle relation when the original relation is unsupported and only one radius
parameter remains. Three-angle cosine Gram consequences with total angle
multiplicity at most six are then constructed in that one parameter. The
original multivariable expansion bounds remain in place. These are necessary
equations only; every resulting candidate still requires exact verification of
all original faces, vertex fans, and convexity conditions.

This resolves body 153 (J48) through `--radius-only`. Its postponed belt fan
has angle relation `3*A + B + C = 2*pi`. Substitution reduces its closure to
`5*d^4 + (15+7*sqrt(5))*d^2 - 9 - 4*sqrt(5) = 0`, where `d=s1-s0`.
The common squared-radius field has degree eight. The implementation contains
no body identifier, precomputed root, or stored answer polynomial.

## Exact roots with algebraic coefficients

Version 3.0.2 uses the proper field norm when polynomial coefficients contain
algebraic constants such as sqrt(5). The resulting rational polynomial is
factored and its real roots isolated. Each candidate must satisfy the original
polynomial in the original coefficient embedding exactly, since the norm also
contains roots belonging to conjugate coefficient embeddings. Rational
coefficient polynomials retain their existing extraction path.

This repairs J20 and applies throughout the solver, without body-name tests,
coefficient guessing, or a dependency upgrade.

## Exact constant reduction

Version 3.0.4 fixes two failures after algebraic substitution. Constant
equations are cancelled before testing for zero; unresolved algebraic
identities are checked in their exact number field, including isolated real
roots. A tiny nonzero constant remains a contradiction. Unsupported identities
remain unresolved, and the existing time limits still apply.

Constant numerators and denominators of radicands are handled as coefficients.
Nonconstant radicands are factored using their explicit polynomial variables,
so algebraic coefficients are not mistaken for missing or inferred variables.

These fixes allow body 67's already-derived values 3/4, 7/4, and
(5+2*sqrt(5))/4 to reach independent verification. They apply to the algebraic
operations themselves and contain no body identifiers or stored answers.

## Coplanar face subdivisions

Version 3.0.3 accepts coplanar regular polygons when they exactly tile a full
convex supporting facet. The independent geometric checker groups faces by
their exact supporting plane, orients their boundaries consistently, and
checks that internal edges cancel to one strictly convex boundary containing
all vertices on the plane. This certifies complete coverage without overlapping
face interiors. Nonconvex, overlapping, incomplete, and entirely flat inputs
remain rejected. No floating-point tolerance or body identifier is used.

This fixes census 37, whose center is on a square face and whose boundary
contains two rhombi subdivided into equilateral triangles. Its exact squared
radii are 1/2, 3/2, and 11/6. Both the solver and standalone saved-certificate
checker use the same independent geometric checks. Reports distinguish the
listed regular faces from their merged maximal supporting facets.

## Boundary and exterior sphere centers

Version 3.0.1 retains the supplied or discovered sphere partition when its
center is not strictly interior. It skips unsigned radial equations in that
case and uses exact metric geometry instead. It does not relax the radial
positivity threshold or silently replace the center by the centroid.

This fixes J92: the exact center is the center of its hexagonal face, and the
squared radii are 1 and (3+sqrt(5))/2. The geometric verifier derives its center
classification from exact face-plane signs, independently of numerical input.
Version 3.0.3 also verifies coplanar facet subdivisions through the geometric
route, as described above. The existing interior radial route remains supported.

## Results and certificates

Every run writes `result.json`, `progress.json`, and `run.log`. Completed stages
also write normalized input with source hashes, exact equations, and derivation
traces (`local-derivation.json`, `gram-derivation.json`, and `derivation.json`). The radius fallback records
its angle equations in `equations.json`; the main path records its Gram equations.
The result records whether `radius_only` was requested. A successful full
radius-only calculation reports `search_route: radial_only`; a successful
short prepass retains `search_route: local_radius`.

A successful run additionally writes `radii.json`, `certificate.json`, and
`verification.json`. The geometric route includes the exact coordinate metric,
vertex coordinates, and sphere center in `radii.json`, so the standalone checker
can independently reconstruct every geometric obligation. Exact squared radii are rational functions of a real root
`alpha` of an irreducible polynomial, selected by a rational isolating interval.
The requested polynomial is

\[
Q(T)=\prod_j(T-s_j),\qquad s_j=r_j^2.
\]

`Q_coefficients` gives its coefficients in descending powers of T. `radii` gives
the positive square roots. The degree of the primitive field polynomial need
not equal the number of sphere groups. For rational data the primitive element
is harmlessly represented by `alpha=0`.

`--time-limit` is a hard solving wall-clock limit, including verification. A
timeout or interruption gives `UNRESOLVED`. `--overwrite` replaces this run's
previous output files; without it an existing result is protected.
`--equations-only` saves exact radial equations and deferred closure templates
without attempting to solve them.

## Short field expressions (v3.0.8)

When a closure parameter already generates the field containing the radius
coefficients, the solver retains it. For example, a coefficient such as
sqrt(5) may belong to Q(d) even though it is initially written separately.
An exact coefficient-field embedding and a rational-span closure check prove
that the radii generate precisely that field within a two-second optional
budget; otherwise the previous primitive-element construction is used.
The existing rational-coefficient fast path is unchanged.

After certification, a separate presentation stage writes:

- `presentation.html`: an offline report with readable equations.
- `presentation.md`: the same equations in Markdown/LaTeX.
- `presentation.json`: the generator, shared quantities, short rational
  formulas, exact conversion map, and equivalence-check results.
- `presentation.log`: completed improvements.

The stage reconstructs rational functions modulo the defining polynomial,
tries generators drawn from radius differences and individual radii, and
factors out a common squared-radius offset. It also shares affine relations
between the coefficients in expressions of the form a+b*t+c. Only a shorter
complete presentation replaces the previous one. No body labels, stored
answer polynomials, or recognition from numerical approximations are used.

Every replacement is checked exactly against the certified field, including
nonzero denominators and the real embedding. This is a bounded readability
heuristic, not a guarantee of the shortest possible formulas. It incorporates
the earlier v4 presentation experiment while retaining the current v3 solver.
The proof representation stays in `radii.json`; Q is displayed as a product
of the already defined squared-radius factors.

The default presentation search budget is five seconds, separate from solving.
Report rendering follows that search. Use `--presentation-time-limit 10` to
increase it, or `--no-presentation` to omit the stage. Both options are accepted
by `run_body.py`, `solve.py`, `run.sh`, and `benchmark.py`. A search timeout or
unsupported simplification preserves the completed SOLVED certificate.

Existing certified results can be reformatted without rerunning the solver:

```sh
/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 present.py \
  /path/to/existing/certified/output --out out/readable --time-limit 5
```

This writes only to the chosen destination. To check the change independently:

```sh
/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 presentation_check.py \
  /path/to/existing/certified/output --presentation out/readable/presentation.json
```

`check.py` additionally checks `presentation.json` when it accompanies the
certificate. A tampered display fails that combined check.

## Independent verification and tests

```sh
./run.sh --help
# Use the same Python selected by run.sh for these two commands:
/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 check.py out/J89
/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 benchmark.py \
  --out out/worked-examples --time-limit 120 --overwrite --recheck
```

The standalone checker uses the original faces, vertex groups, and saved exact
algebraic data. For geometric certificates these include the coordinates,
positive metric, and center. It does not read BI coordinates as evidence or
trust the solver's elimination trace. Certificate expressions use a restricted arithmetic parser.

See `VALIDATION.md` for measured results from this version. Earlier solver
versions are separate packages and are unchanged.
