# Admissible full belts around a regular heptagon

## Result and assumptions

Exactly two admissible full belts occur, up to congruence and reversal:

\[
\boxed{((7,4,4))^7},\qquad
\boxed{((7,3,3,3))^7}.
\]

They are the seven-square prism belt and the fourteen-triangle antiprism
belt. Both have fourteen vertices, all on one sphere. Both have an actual
convex completion by the opposite regular unit heptagon.

The distinguished heptagon is included in the geometric fragment. Faces
are prescribed regular unit-edge polygons with proper edge-to-edge
intersections: two distinct faces intersect in at most one common vertex
or one common edge. Every prescribed vertex of a completion has strictly
positive angular defect. Adjacent faces may be coplanar. A full belt
includes every face incident with a base vertex, including faces meeting
the base only at that vertex.

Admissibility requires mutually compatible regular-face completions of
the open stars, supporting the **whole augmented fragment**. The
exclusions below use necessary conditions for such a finite compatible
completion; the survivors are exhibited as complete convex bodies. No
completed-polyhedron catalogue is used. No largest-face assumption or
artificial bound on future face sizes is used. In particular, the outer
boundary is not assumed to be a simple ring, nor its vertices distinct.

## 1. Exact local alphabet

Write \(\alpha_j=\pi-2\pi/j\). At a base vertex the sum \(A\) of the
other face angles satisfies

\[
5\pi/7<A<9\pi/7.
\]

The first inequality is the strict angle polygon inequality for a
nondegenerate convex vertex; the second is positive defect. Four other
faces already contribute at least \(4\pi/3>9\pi/7\). Exact rational
inequalities therefore give the following exhaustive marked links,
modulo reversal:

- \((7,3,k)\), \(4\le k\le41\);
- \((7,4,k)\), \(4\le k\le9\);
- \((7,5,5)\), \((7,5,6)\);
- \(A=(7,3,3,3)\);
- \(Q_1=(7,3,4,3)\), \(Q_2=(7,3,3,4)\);
- \(P_1=(7,3,5,3)\), \(P_2=(7,3,3,5)\).

There are 51 links modulo reversal and 97 directed links. These are
local possibilities, not admissible belts. Their enumeration and the
supported transition equations are in [local-transitions.md](local-transitions.md).

In each four-face link, the middle face of the three displayed after
the heptagon is vertex-only. Thus \(Q_1\) has a vertex-only square and
\(P_1\) a vertex-only pentagon. The two \(Q_2\) or \(P_2\) directions
are retained when gluing.

## 2. Geometric states and the mixed-family exclusions

Put \(\theta=\pi/7\), \(C=\cos(2\theta)\), \(S=\sin(2\theta)\),
\(R=1/(2\sin\theta)\), and \(m=\cot\theta/2\). A triangle on a base
edge has vertices

\[
(m,-1/2,0),\quad(m,1/2,0),\quad(a,0,z),
\qquad z^2=3/4-(a-m)^2.
\]

Its exposed angle \(\beta\) has \(\cos\beta=1/2-Sa\). Reflection
in the base edge's perpendicular-bisector plane preserves both the base
and the attached regular face. Hence its exposed angle is the same at
both endpoints. This gives edge consistency without any assumption on
the remote boundary.

The Gram equations, together with support against the retained base,
select a unique branch at each four-face link. Links \(A,Q_1,P_1\)
give involutions \(f_3,f_4,f_5\) on the triangle state, with unique fixed
radii \(R,\sqrt2R,\varphi R\), respectively, where
\(\varphi=(1+\sqrt5)/2\). The \(Q_2,P_2\) transitions are unique in
either direction between triangle and square/pentagon states. Coplanar
endpoints remain in the supported intervals.

Unlike the nonagon, the heptagon has overlapping state intervals: \(A\)
can potentially neighbor \(Q_1\), and \(P_1\) can potentially neighbor
\(Q_1,Q_2\), across an edge triangle. These interfaces require separate
exclusions before using odd-cycle rigidity.

The following exact checks provide them:

| Interface | Necessary completion and obstruction |
|---|---|
| \(P_1\)–trivalent \((7,3,k)\) | The angle budget permits no new face. Closure would require the same triangle/\(k\)-gon dihedral with third angle \(\alpha_7\) at one end and \(\alpha_5\) at the other, which is impossible. |
| \(P_1\)–\(Q_1\) | Only one triangle or square fits the remaining budget, but the gap is strictly greater than \(90^\circ\). |
| \(P_1\)–\(Q_2\) | A heptagon is forced. An exhaustive interval cover encloses every possible fold; the neighboring triangle then requires two triangles, and both possible completions have a face plane cutting the existing fragment. |
| \(A\)–\(Q_1\) | The other triangle at \(Q_1\) forces a \(Q_1\)–\(Q_2\) pair with a gap strictly between \(60^\circ\) and \(90^\circ\). Only two triangles can fill it; both roots have a required face plane cutting the retained base or an existing square corner. |

The complete arguments, support witnesses and exact interval certificates
are in [pentagon-proof.md](pentagon-proof.md) and
[links7-bridge.md](links7-bridge.md). The latter also uses the finite
trivalent endpoint exclusions in [quad-proof.md](quad-proof.md): a
\(Q_1\) or \(Q_2\) triangle can meet a trivalent \((7,3,k)\) endpoint
only at \(k=7\).

State intervals rule out the remaining direct interfaces. Reflection
also prevents a trivalent intervening face from switching its neighboring
face size between its two base endpoints. The three families may now be
treated separately.

## 3. Elementary components and the all-triangle family

The link \((7,4,4)\) propagates uniquely: the exposed square angle
\(90^\circ\) lies below the \(Q_2\) square range. Its common unit side
edges are perpendicular to both neighboring base edges, so the geometry
is exactly the square prism belt.

A component containing only trivalent links satisfies
\(s_{i-1}=s_{i+1}\) for its cyclic base-edge face sizes, by reflection.
On the odd seven-cycle this forces all sizes equal. The local alphabet
allows only all squares or all pentagons. The all-pentagon case is
impossible: each outer pentagon–pentagon shoulder forces a heptagon, and
two such distinct heptagons meet the central outer vertex of a lateral
pentagon, contributing

\[
\alpha_5+2\alpha_7=71\pi/35>2\pi.
\]

Proper face intersections prevent the two forced heptagons from being
the same face. This also deals with possible remote vertex identifications.

After the interfaces above are removed, the \(A\) family has edge
triangles and isolated pentagons or hexagons. Each nontriangle block
uses two base vertices and preserves the adjoining triangle state.
The remaining number of \(A\) transitions is odd, so closure under
\(f_3\) forces its fixed radius \(a=R\).

A pentagon would instead require \(Sa=\cos(\pi/5)\), and a hexagon
would require \(Sa=1\). Neither equals \(SR=\cos(\pi/7)\), which lies
strictly between them. Only the constant-radius all-triangle belt
remains: the regular heptagonal antiprism belt.

## 4. Excluding the square-connector family

The remaining links in a component containing \(Q_1,Q_2\) are these
two links, trivalent triangle–heptagon links, and trivalent square–pentagon
or square–hexagon links. Without triangles the edge labels alternate,
impossible on seven edges.

Every maximal nontriangle block is a single heptagon or
\(SXSX\cdots S\), where \(X\) denotes pentagons or hexagons. It
consumes an even number of base vertices and preserves the triangle
state. The remaining \(Q_1\) count is odd, so \(f_4\) forces
\(a=\sqrt2R\). Exact state comparisons exclude the heptagon and
pentagon/hexagon blocks at this fold. Only \(Q_1,Q_2\) remain, with
at least one \(Q_1\).

At the fixed fold a \(Q_1\)–\(Q_1\) triangle apex has a gap strictly
between \(90^\circ\) and \(108^\circ\); none of the possible single
triangle, square or pentagon completions fits. Every \(Q_1\) therefore
has \(Q_2\) across its two edge triangles. Each mixed apex forces a
heptagon.

At the opposite corner of the \(Q_1\) square, the square and these two
distinct heptagons leave only \(90/7\) degrees, less than a triangle
angle. Their remaining edge rays must therefore coincide. They do not:
with \(p=\sqrt{1+C}\), their transverse components are \(Y,-Y\), where

\[
Y=-\frac{\sqrt2(2p^3-2p-1)}{2(2p+1)}\ne0.
\]

Indeed \(C>9/16\) gives \(p>5/4\), and the numerator polynomial is
positive there. Thus the two forced faces are incompatible. Full exact
derivations and symbolic verification appear in [quad-proof.md](quad-proof.md).

## 5. Excluding the pentagon-connector family

After the mixed interfaces are removed, nontriangle blocks are a single
\(K\) or \(P(HP)^j\), where \(P\) is a pentagon and \(H\) a hexagon.
Each consumes an even number of base vertices and preserves triangle
state. The remaining \(P_1\) count is odd, so \(f_5\) forces
\(a=\varphi R\).

At this fold the triangle exposed angle is strictly between
\(\alpha_{21}\) and \(\alpha_{22}\), excluding every integer \(K\).
The pentagon fold differs from the one forced by a neighboring hexagon.
Only \(P_1,P_2\) remain.

A \(P_1\)–\(P_1\) triangle apex permits only one further triangle,
but its gap exceeds \(60^\circ\). A \(P_2\)–\(P_2\) triangle apex
has gap strictly between \(\alpha_{27}\) and \(\alpha_{28}\); no single
regular polygon fits. Its remaining positive angle budget is less than
\(180^\circ\), so the only possible multiface completion is a triangle
and a pentagon.

One geometric root for that fan cuts the existing fragment. The other
forces a new pentagon whose neighboring old pentagon shoulder then has
existing angles \(108+60+108=276\) degrees. Only a triangle could complete
that shoulder, but its gap is strictly between \(60^\circ\) and
\(90^\circ\). Reflection covers the other face order. Therefore
\(P_2\)–\(P_2\) across a triangle is impossible too.

Every triangle must consequently join \(P_1\) to \(P_2\), while each
edge pentagon joins two \(P_2\) vertices. The cyclic link word must
repeat \(P_1,P_2,P_2\). Since \(3\nmid7\), no such belt closes.
The gap and supporting-plane certificates are detailed in
[pentagon-proof.md](pentagon-proof.md).

## 6. Explicit realizations, spheres and curvature

With indices modulo seven, put

\[
V_i=(R\cos(2\pi i/7),R\sin(2\pi i/7),0).
\]

For the prism take \(U_i=V_i+(0,0,1)\), and the seven squares
\((V_i,V_{i+1},U_{i+1},U_i)\).

For the antiprism take

\[
U_i=(R\cos((2i+1)\pi/7),R\sin((2i+1)\pi/7),H),
\qquad H^2=1-\frac1{4\cos^2(\pi/14)},
\]

and the fourteen triangles
\((V_i,V_{i+1},U_i)\), \((V_i,U_{i-1},U_i)\).
These include the seven vertex-only triangles required by the definition
of a full belt. All edges have length one. Adding the opposite heptagon
on the \(U_i\) closes a convex prescribed surface in either construction.

All fourteen vertices lie on

\[
x^2+y^2+(z-H/2)^2=R^2+H^2/4,
\]

using \(H=1\) for the prism. Each vertex of either completed surface has
defect \(2\pi/7\), so the total is exactly \(4\pi\). Equality is allowed
by the curvature bound “does not exceed \(4\pi\).”

## Verification

The enumeration uses exact fractions. The Gram residuals and fixed-fold
identities have symbolic checks. Inequality certificates use rational
outward intervals, including whole-interval covers where folds vary;
the displayed decimal approximations do not decide any classification.
Exploratory files with `probe` in their names are not proof certificates.
The supporting notes record all certificate files. The component logic
and certificate calculations were also checked independently; see
[independent-audit.md](independent-audit.md).
