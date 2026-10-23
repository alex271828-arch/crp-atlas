# Admissible full belts around a regular nonagon

## Result and scope

Exactly two admissible full belts occur, up to congruence and reversal:

\[
\boxed{((9,4,4))^9},\qquad \boxed{((9,3,3,3))^9}.
\]

They are the nine-square prism belt and the eighteen-triangle antiprism belt.
The distinguished regular unit-edge nonagon is included. Each belt has eighteen
vertices, all on one sphere, and extends by adding the opposite regular nonagon.

The assumptions are convexity, regular unit-edge prescribed faces, proper
edge-to-edge face intersections, and strictly positive defect at prescribed
vertices. Coplanar neighboring faces are allowed. An admissible completion of
the original open vertices must use mutually compatible actual polygons whose
planes support the whole augmented fragment. The proof uses necessary conditions
for such a finite completion; it does not assume a completed body or identify
independent convex vertex cones with admissibility. No ordinary outer ring or
distinctness of all outer vertices is assumed in the exclusions.

No complete-polyhedron catalogue, maximum-face hypothesis, or artificial bound
on future face sizes is used.

## 1. Exact local links

Set alpha_k = pi - 2pi/k. At a base vertex the sum A of the other angles satisfies

\[
7\pi/9<A<11\pi/9.
\]

The complete list, with the base marked and reversal identified, is

* (9,3,k), 4 <= k <= 17;
* (9,4,k), 4 <= k <= 7;
* (9,5,5);
* (9,3,3,3);
* Q1 = (9,3,4,3);
* Q2 = (9,3,3,4), with its reverse representing the other direction.

There are 22 links modulo reversal and 40 directed links. These are angle
possibilities, not yet belts. The exact reciprocal inequalities and derived
finite bounds are checked by [links9-enumerate.py](links9-enumerate.py).

In Q1 both base-edge faces are triangles and a square meets the base only at the
vertex. In Q2 a square is on a base edge and a triangle meets the base only at
the vertex. In the all-triangle link a third triangle meets the base only at
the vertex.

## 2. Supported states and edge consistency

Put theta = pi/9, C = cos(2theta), S = sin(2theta),
R = 1/(2sin theta), m = 1/(2tan theta).
For a base-edge triangle its apex has radius and height

\[
a=m+(\sqrt3/2)\cos t,\qquad u=(\sqrt3/2)\sin t.
\]

The exposed angle beta between the unused base ray and the outgoing triangle
ray satisfies cos beta = 1/2 - Sa. Reflection in the perpendicular-bisector
plane of a base edge preserves the base and its regular edge face, so the
exposed angle and the edge-face fold agree at both endpoints.

The supporting-plane calculations from the general triangular lemma and the
Q1/Q2 Gram identities apply at theta = pi/9. They retain coplanar endpoints and
give:

| Link and edge face | Exposed-angle interval | Single polygons in that interval |
|---|---|---|
| (333), triangle | [acos(1/2-C), 120 degrees] | 5,6 |
| Q1, triangle | [acos(sqrt(3)/2-2C), 150 degrees] | 8,9,10,11,12 |
| Q2, triangle | [140 degrees,150 degrees] | 9,10,11,12 |
| Q2, square | [acos(1-sqrt(3)C),120 degrees] | 6 |

The triangle intervals for (333) and (334) are disjoint. Exact rational
enclosures for all separating comparisons are checked by
[links9-state-certificate.py](links9-state-certificate.py). The general Gram
identities are in [the earlier local derivation](../belt-study-11-plus/n11-proof.md#7-explicit-q2-gram-certificate-for-the-local-intervals)
and [the square-connector derivation](../belt-study-11-plus/n11-square-connector.md#1-exact-range-imposed-by-the-supporting-planes).
These identities are polynomial in C and remain valid at the nonagon value;
the required signs hold throughout the displayed intervals.

At (333), the supported transition from one edge triangle to the other is a
unique involution f3, with fixed radius R. At Q1 the supported transition is a
unique involution f4, with fixed radius sqrt(2)R. For f4, equal-radius adjacent
apices have distance sqrt(2), the diagonal of the connector square; this gives
the fixed radius. Q2 gives a unique supported transition in either direction
between a triangle fold and a square fold. The support signs select the unique
circle-intersection branch. These are necessary local conditions from support
of the retained base; requiring more supporting inequalities can only discard
candidates.

## 3. The elementary components

The link (9,4,4) cannot transition to Q2 across a square: its exposed angle is
90 degrees, below the Q2 range. It propagates to the square prism belt. At every
base vertex the shared square side edge is perpendicular to both base edges,
so all these unit side edges are perpendicular to the base. This proves its
unique geometry.

The links (9,4,5) and (9,4,7) force pure alternation around the base, impossible
on nine edges. They cannot reach Q2 because 108 degrees and alpha_7 lie outside
the Q2 square interval. Likewise (9,3,k), k = 4,7,13,...,17, force pure
triangle/k-gon alternation and are impossible on the odd base cycle.

The link (9,5,5) propagates to nine pentagon edge faces. At the outer endpoint
of a shared pentagon-pentagon edge, the existing angle sum is 216 degrees and
the free-ray angle is 140 degrees, by reflection along that shared edge. The
remaining angle sum must be below 144 degrees. Two triangles total only 120
degrees; every other multi-face collection exceeds that budget. Hence exactly
one nonagon is required at each such shoulder.

For a pentagon with cyclic vertices A,B,U,X,V and base edge AB, the shoulders
U,V therefore force completing nonagons along UX and VX. They are distinct:
one prescribed face cannot share two edges with the pentagon. At X their angles
and the pentagon angle total 140+140+108 = 388 degrees, impossible. Thus the
all-pentagon component is excluded by compatible completion.

## 4. The (333) component

Its base-edge labels are triangles, with isolated pentagons or hexagons. Each
pentagon or hexagon has trivalent endpoints and uses two base vertices; crossing
this block preserves the neighboring triangle fold. All remaining base vertices
apply f3. Their number is 9 minus twice the number of these nontriangle edges,
and hence is odd. Closing the cycle forces the fixed state a=R everywhere.

A pentagon endpoint instead requires a=D/S, where D=(1+sqrt(5))/4=cos(pi/5).
This is less than R, since SR=cos(pi/9)>D. A hexagon requires a=1/S>R.
Consequently neither pentagons nor hexagons can occur.

The remaining all-triangle belt has constant radius R and constant positive
height, so it is the regular nonagonal antiprism belt. The odd-cycle argument
already fixes its geometry; no separate even-cycle alternating realization is
left at n=9.

## 5. The (334) component and odd-cycle rigidity

The remaining base-edge labels are T=triangle, S=square, H=hexagon, and K of
size 8,...,12. Its links are Q1, Q2, (9,4,6), and (9,3,k) with 8<=k<=12.
If T is absent the word alternates S,H, impossible on an odd cycle. Otherwise
each maximal non-T block is either K or S(HS)^j, j>=0.

Every block consumes an even number of base vertices and preserves the
neighboring triangle fold. The other vertices apply f4. Their number is odd,
so closing the nine-cycle forces the f4 fixed state

\[
a=\sqrt2 R
\]

on every triangle edge. In particular at least one Q1 occurs.

A K endpoint of size k requires

\[
a_k=\frac{1/2+\cos(2\pi/k)}{S}.
\]

The exact state certificate proves a_10 < sqrt(2)R < a_11. Monotonicity in k
therefore excludes every integer K size in the range 8,...,12. A block containing
H forces the coplanar-two-triangle endpoint of Q2, hence radius a_9, also
different from sqrt(2)R. Thus only Q1 and Q2 remain, at the fixed triangle fold.

## 6. Compatible completion excludes the remaining Q1/Q2 patterns

Write p=sqrt(1+C). When both endpoints of a base-edge triangle are Q1, the
outer apex has one triangle and two squares, total 240 degrees. The gap gamma
satisfies

\[
\sin(\gamma/2)=p-1/2>\cos(\pi/5)=\sin54^\circ.
\]

The strict inequality is the same certified comparison a_10<sqrt(2)R used
above. Hence gamma>108 degrees. Under the strict remaining budget of 120
degrees, only one triangle, square, or pentagon could fit, and none spans the
gap. Therefore a Q1 cannot neighbor another Q1 across a triangle edge.

Every Q1 consequently has Q2 at both neighboring base vertices. At a triangle
with Q1 at one endpoint and Q2 at the other, the existing faces at the outer
apex are triangle, square, triangle, total 210 degrees. The exact supported
geometry gives

\[
\cos\gamma=-C,\qquad\gamma=140^\circ.
\]

The remaining angle sum must be strictly below 150 degrees. Two triangles
cannot span the gap; a triangle and a square would give zero defect. Hence a
single nonagon is forced.

To verify the geometric identity directly, align the base-edge triangle as

\[
A=(m,-1/2,0),\quad B=(m,1/2,0),\quad T=(r,0,z),
\]

where r=p/S, m=p^2/(2S), x=r-m and z^2=3/4-x^2. The supported Q2 square fold is

\[
q=\frac{p^2}{S(2p+1)},\qquad
w=\frac{2(p+1)z}{2p+1}.
\]

The outgoing unit rays from T belonging to the Q1 square and the Q2 triangle are

\[
v=(rC-m,\,1/2-p,\,z),\quad
u=(Cq-x,\,1/2+Sq,\,w-z).
\]

Substitution of C=p^2-1 and S^2=p^2(2-p^2) gives v.u=-C. The support signs
selecting the displayed Q2 fold are positive, so the other algebraic branch
does not evade this conclusion. Full derivation and symbolic checks are in
[quad-proof.md](quad-proof.md).

Now take the vertex-only square at any Q1. At each of its two corners adjacent
to the base vertex, the preceding mixed-apex argument forces a completing
nonagon. Both nonagons meet its opposite corner, and they are distinct because
one face cannot share two square edges. At that corner the angles total

\[
90^\circ+140^\circ+140^\circ=370^\circ,
\]

which is impossible. Thus no Q1 can occur, contradicting its required odd
number. The entire (334) component is excluded. This argument also handles
the apparently possible three-periodic Q1,Q2,Q2 word; divisibility of nine by
three does not make its outer completions compatible.

## 7. Explicit survivors, curvature, and sphere

With R=1/(2sin(pi/9)), let

\[
V_i=(R\cos(2\pi i/9),R\sin(2\pi i/9),0),\quad i\pmod9.
\]

For the prism set U_i=V_i+(0,0,1), and take the nine squares
(V_i,V_{i+1},U_{i+1},U_i).

For the antiprism set

\[
U_i=(R\cos((2i+1)\pi/9),R\sin((2i+1)\pi/9),H),
\qquad H^2=1-\frac1{4\cos^2(\pi/18)},
\]

and take triangles (V_i,V_{i+1},U_i), (V_i,U_{i-1},U_i). The U vertices form a
regular unit-edge nonagon in both cases. Adding it closes a convex prescribed
surface; each of its eighteen vertices has defect 2pi/9. The total curvature
is exactly 4pi, as allowed by the bound "does not exceed 4pi."

All belt vertices lie on the sphere

\[
x^2+y^2+(z-H/2)^2=R^2+H^2/4,
\]

where H=1 for the prism. The coordinate construction proves all eighteen
vertices are distinct. Completeness of the exclusions did not assume this
in advance.

## Verification status

The link enumeration and state separation use exact rational arithmetic and
outward rational enclosures of algebraic/trigonometric quantities. The mixed
Q1/Q2 geometry is checked by symbolic polynomial identities with positive
support signs. Numerical exploration is not used as a completeness argument.
The classification-free primary face-restriction paper was consulted for
context, but no catalogue or black-box nonagon classification is used here.
