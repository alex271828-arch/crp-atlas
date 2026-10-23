# Independent review of the n>=11 belt argument

Reviewed the final n=11 state proof, its two local completion exclusions, and
the extension of the triangular rigidity argument to every n>=12. No
substantive gap remains in these arguments under the stated proper-intersection,
whole-fragment convexity, and joint regular-completion assumptions.

## n=11 local states and the odd-cycle argument

The exact enumeration gives 16 marked links modulo reversal and 29 directed
links. In particular it retains (11,3,4), (11,4,6), and both marked orders of
the (3,3,4) sectors, as required; the n>=12 reductions cannot simply be applied
at n=11.

The two Q2 Gram calculations in Section 7 of `n11-proof.md` were independently
checked. Their support-selected circle intersections are unique in **both**
directions, triangle to square and square to triangle. The factorizations

\[
-3C(C+k)(1-k^2)/4,
\qquad
(1-k^2)((1-k)^2-3C^2)/4
\]

are correct. The intervals preserve equality, including the conditional-edge
states. The Q1 Gram calculation in `n11-square-connector.md` is also correct.
These unique symmetric transitions give involutions; the parity argument does
not require an unproved global monotonic-branch assertion.

The exposed-angle intervals separate the elementary, (333), and (334)
components exactly as claimed. Their strict comparisons have been checked by
rational interval arithmetic. Every non-triangular block in the (334)
component is either K or S(HS)^r, and consumes an even number of base vertices.
Reflection and transition uniqueness imply equal triangle states on its two
sides. Therefore an odd 11-cycle has an odd number of Q1 involutions and
forces their fixed state r=sqrt(2)R. This counts base vertices, not presumed
outer-ring vertices.

The fixed state excludes K of either size 11 or 12. Indeed, after multiplying
by S the three radii are 1/2+C, sqrt(1+C), and (1+sqrt(3))/2; their strict
ordering follows from C<sqrt(3)/2. A hexagon block forces the coplanar Q2
endpoint and the first of these radii, and is therefore excluded too. The
(333)/hexagon component has the same valid odd-involution argument, with
fixed radius R instead of 1/S.

## The two remaining local completion exclusions

For Q1--triangle--Q1, the Gram interval forces a free-ray gap greater than
108 degrees. The three existing sectors have total 240 degrees, leaving
strictly fewer than 120 degrees. Only one triangle, square, or pentagon could
fit, and none spans that gap. The proof uses necessary supporting-plane tests
against the existing base, and remains valid with coplanar edges.

For Q2--triangle--Q2 at the fixed state, I independently rederived the local
coordinates and the two supporting tests in `n11-q2-pair.md`. With
L=Cx-S/2, those tests are exactly

\[
qz-Lw,\qquad wx-Cqz.
\]

The discarded circle-intersection solution has opposite signs, hence cuts
the existing base. The retained solution gives d=1/2+Sq and the rigorous
certificate encloses it strictly between cos(pi/34) and cos(pi/35). Thus the
missing angle lies strictly between alpha_34 and alpha_35. No single regular
polygon of any integer size fits; two or more faces have total at most 168
degrees under the strict 180-degree budget, and also cannot fit.

The residual link cycle must consequently repeat Q1,Q2,Q2: every square edge
joins two Q2 vertices, while every triangle edge joins unlike Q-types. Since
3 does not divide 11, this last component is impossible. Possible additional
outer incidences cannot repair either local gap without obeying the same
face-angle budget and supporting-plane constraints.

## Uniform triangular rigidity for n>=12

The generalization in `general-tri-proof.md` is valid with strict endpoint
inequalities. At the higher parity the gap exceeds alpha_n>=150 degrees, so
the only multiple-face option is (3,5). At the lower parity its gap exceeds
pi-4pi/n>=120 degrees, and consistency with its neighbors again forces (3,5).
The exact absolute-height identity then puts a required completing triangle
apex below the base: its bracket is less than

\[
(D-\cos^2(\pi/n))/\sin(2\pi/n)<0,
\qquad D=(1+\sqrt5)/4.
\]

This uses whole-fragment supporting planes and agreement of neighboring
faces, not independent convex cones. It treats n=12's lower-bound equality
correctly, because the actual gap is strictly larger than that bound.

## Executed checks

The following exact checks passed during this review:

- `verify_links11.py`;
- `verify_n11_gram.py` (polynomial identities by exact interpolation and
  explicit degree bounds);
- `verify-n11-state-separation.py` (independent interval checks of the state
  separation and radius ordering);
- `certify_n11_q2_pair.py` (rational interval enclosure and support signs);
- `verify-general-local.py` (uniform n>=12 alphabet and hexagon angle gaps);
- `../belt-study-12/verify_triangular_identities.py` (symbolic identities with
  arbitrary trigonometric parameter, including the completion-height equation).

The symbolic script uses the existing local SymPy dependency through
`PYTHONPATH=outputs/.research-deps`. Numerical probes are not used as proof.
The retained n=12 square--pentagon exclusion is the existing executed exact
certificate from the dodecagonal study; this review does not replace that
certificate.
