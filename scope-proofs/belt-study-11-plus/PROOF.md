# All admissible belts around a regular n-gon for n >= 11

## The result, in marked-link notation

For every integer n >= 11, precisely two admissible full belts occur:

\[
\boxed{\bigl((n,4,4)\bigr)^n}
\qquad\text{and}\qquad
\boxed{\bigl((n,3,3,3)\bigr)^n}.
\]

The first is the regular n-square prism belt. The second is the regular 2n-triangle antiprism belt. These statements classify the actual geometric fragments up to congruence, not just their link multisets. Both explicitly extend by an opposite regular n-gon. All their vertices, including the distinguished n-gon, lie on one sphere.

The proof below uses no completed-polyhedron catalogue or restricted census. In particular, no global maximum face size is assumed. The n=11 step permits arbitrarily large future faces and includes coplanar prescribed faces.

## 1. Assumptions and admissibility

The distinguished regular unit-edge n-gon N is retained. Its full belt contains every prescribed face meeting any vertex of N, including vertex-only faces. All N-vertex links are closed. Faces are regular unit polygons with proper edge-to-edge intersections. The fragment and any proposed completion lie on the boundary of a common three-dimensional convex set; every prescribed face plane supports the whole augmented fragment.

At all originally open vertices, a joint completion witness must consist of actual regular polygons closing their links with strictly positive defect. Neighboring completions agree about shared faces. New vertices may remain open: a complete body is not assumed. The exclusions below use necessary conditions for this finite witness. The explicit caps establish actual complete-body extendibility for the two survivors.

No annular outer boundary or distinctness of outer vertices is assumed during exclusions. In the surviving square belt, distinctness follows from the perpendicular unit side edges. In the surviving triangular belt, a coincident pair of base-edge apices would be unit distance from three base vertices, impossible because their circumradius R=(2sin(pi/n))^{-1}>1.

For curvature, use the bound **does not exceed 4pi**. A witness gives the sum of the final defects at the original fragment vertices; a full extension has total4pi. Infima over unrestricted compatible witnesses need not be attained, so they must not automatically be called minima. The proof here needs no positive minimum defect for an unrestricted face alphabet. Both displayed cap constructions have exactly2n vertices of defect2pi/n and attain total4pi.

## 2. The uniform reduction for n >= 12

Let alpha_k=pi−2pi/k. At a base vertex, the sum S of the other face angles satisfies alpha_n<S<2pi−alpha_n. The complete local list for n>=12 is

* (n,3,k), with 6n/(n+6)<k<6n/(n−6);
* (n,4,4);
* (n,4,5), only when12<=n<=19;
* (n,3,3,3).

These strict inequalities prove the finite local list; they are not computational cutoffs on future faces. The exact table and all uniform exclusions are in [general-local-proof.md](general-local-proof.md).

The reductions are:

1. A triangle with a k-gon sector at one base endpoint has the same exposed-ray angle at the other endpoint, by reflection of N and that edge triangle. For k>=7, two triangles cannot span the angle, so another k-gon is forced. At the triangle apex no additional regular face fits positively. Equality of the dihedral along a shared triangle–k-gon edge forces k=n, which contradicts positive defect for n>=12.
2. For a hexagon, the second endpoint has either another hexagon or two coplanar triangles. The first alternative has the same dihedral contradiction. In the second the required apex angle has cosine1/2+C−2C², C=cos(2pi/n). It is strictly between90° and108° for12<=n<=17, and between108° and120° for n>=18. No permitted completing polygon fits. Exact rational bounds certify the separating comparisons.
3. A pentagon–triangle–pentagon sequence has an apex gap too large for the only possible additional triangle. A pentagon–triangle–triangle transition would force triangle apex radius D/sin(2pi/n), D=(1+sqrt5)/4, below the necessary convex support bound cot(2pi/n), since D<cos(2pi/n) for n>10. This excludes every triangle/pentagon mixture.
4. A square–pentagon link forces alternation around N. Odd n is impossible. For even n>12, each outer shoulder requires a single n-gon, because its gap exceeds150° while its remaining angle budget is less than162°. Two distinct such faces at the middle outer vertex of a pentagon exceed360°. The n=12 equality case is excluded by the previously established [exact rational-interval certificate](../belt-study-12/certify_alternating45.py), including its coplanar triangle/square alternative. See [square-pentagon-proof.md](square-pentagon-proof.md).

Thus only the all-square and all-triangle patterns remain. Squares have perpendicular unit side edges and therefore the unique regular prism geometry.

The triangular geometry needs an additional argument. Its adjacent apex parameters lie on a uniquely determined involution, so the radii/heights are constant or alternate between(a,u) and(b,v). In a nonconstant admissibly completable realization, a<R<b and u>v. At the higher, smaller-radius vertices, any single completing face has a plane cutting the existing belt. The remaining possibilities force triangle+pentagon completions; neighboring compatibility forces the same pair at the other parity. Their required triangle apex W=(X,Y,Z) then satisfies

\[
Z(u-v)=(b-a)\left[\frac{DC}{2S^2}\left(\frac1a+\frac1b\right)-m\right],
\]

where theta=pi/n, C=cos2theta, S=sin2theta, m=(2tan theta)^{-1}, D=(1+sqrt5)/4. The necessary strict radius bounds a,b>C/S make the bracket less than

\[
\frac{D}{S}-m=\frac{D-\cos^2\theta}{S}<0.
\]

Hence Z<0, contradicting N's supporting plane. This argument works for every n>=12; it is not restricted to the earlier n>=42 baseline. The full argument is in [general-tri-proof.md](general-tri-proof.md), using the exact identities and branch proof in the earlier [triangular-fragment note](../belt-study-12/fragment-baseline-notes.md).

## 3. Why n=11 needs a separate proof

There are16 marked links modulo reversal (29 directed links):

\[
(11,3,k)\quad(4\le k\le13),\qquad
(11,4,4),(11,4,5),(11,4,6),(11,3,3,3),
\]

and the two distinct marked orders

\[
Q_1=(11,3,4,3),\qquad Q_2=(11,3,3,4).
\]

In Q1 a square meets N only at the vertex, between two base-edge triangles. In Q2 a square lies on a base edge, with a vertex-only triangle beside it. These extra possibilities are retained explicitly; they cannot be removed using the n>=12 local list.

For an edge face, its fold and its exposed-ray angle are identical at its two base endpoints. The exact supporting-plane ranges in [n11-proof.md](n11-proof.md) divide the link possibilities into the following components:

* all-square links, which force the prism;
* alternating3/k or4/5 components, which fail because11 is odd;
* links(11,3,3,3) mixed with(11,3,6);
* a remaining component formed from Q1, Q2, (11,4,6), (11,3,11), and(11,3,12).

The interval derivations use exact Gram determinants and supporting signs, not numerical sampling. They include all coplanar endpoints. In particular they prove that no other transitions between the components are geometrically possible.

## 4. Odd-cycle propagation

For the (333)/hexagon component, crossing a hexagon edge and its two endpoints preserves the adjacent triangle fold. A (333) base vertex applies a unique involution f3 to the fold; its fixed point is the regular antiprism radius R. There are an odd number of these single-vertex transitions on an11-cycle, since each hexagon block uses two vertices. Closure therefore forces the fixed point. A hexagon instead requires apex radius csc(2pi/11)>R, so no hexagon occurs. The regular antiprism remains.

For the other component, label base-edge faces T=triangle, S=square, H=hexagon, K=11-gon or12-gon. If T is absent, S and H would alternate, impossible on an11-cycle. Between consecutive T edges, every non-T block is either

\[
K\qquad\text{or}\qquad S(HS)^j\quad(j\ge0).
\]

Each block uses an even number of base vertices and preserves the triangle fold. Every remaining Q1 vertex applies one common involution f4, whose fixed point has triangle-apex radius

\[
r=\sqrt2R.
\]

There must be an odd number of Q1 transitions. Consequently closure forces this same fixed state on every T edge.

A K block requires radius

\[
r_k=\frac{1/2+\cos(2\pi/k)}{\sin(2\pi/11)},
\]

and exactly

\[
r_{11}<\sqrt2R<r_{12}.
\]

Both strict inequalities follow by squaring from cos(2pi/11)<sqrt3/2. Hence no K block occurs. A block containing H forces the Q2 coplanar-two-triangle endpoint, which also requires r11; it is excluded too. This step explicitly permits, then eliminates, the conditional-edge alternative associated with(11,4,6).

Only Q1 and Q2 remain, in the fixed geometry r=sqrt2R.

## 5. The last Q1/Q2 possibilities

Two Q1 links cannot meet across a base-edge triangle. At its apex, the existing triangle and two squares contribute240°, while the required remaining gap exceeds108°. Positive defect permits only one triangle, square, or pentagon; none fits. The inequality holds throughout the whole convex parameter interval, with an exact supporting-plane proof in [n11-square-connector.md](n11-square-connector.md).

Two Q2 links cannot meet across a base-edge triangle either. At the parity-forced radius r=sqrt2R, its apex already has three triangles. The exact required remaining angle γ satisfies

\[
\alpha_{34}<\gamma<\alpha_{35}.
\]

A single regular polygon of any integer size cannot fit. Multiple completing faces have total below180° and at most168°, also too small. The coordinate derivation is in [n11-q2-pair.md](n11-q2-pair.md). The rational-interval certificate [certify_n11_q2_pair.py](certify_n11_q2_pair.py) verifies the correct convex root and these strict comparisons with outward bounds; all logical decisions use exact rational endpoints.

Thus every Q1 must have Q2 on both sides. Every Q2 must have another Q2 across its square edge, and Q1 across its triangle edge. The link sequence is forced to repeat

\[
Q_1,Q_2,Q_2.
\]

It cannot close after11 vertices. This excludes the last component, completing the n=11 proof without a finite future-face alphabet.

## 6. Explicit geometry, extension, and common sphere

Let the base vertices be

\[
V_i=(R\cos(2\pi i/n),R\sin(2\pi i/n),0),\quad
R=\frac1{2\sin(\pi/n)}.
\]

**Prism belt.** Set U_i=V_i+(0,0,1). Its n faces are squares(V_i,V_{i+1},U_{i+1},U_i). Every base link is(n,4,4).

**Antiprism belt.** Set

\[
U_i=(R\cos((2i+1)\pi/n),R\sin((2i+1)\pi/n),H),
\qquad H^2=1-\frac1{4\cos^2(\pi/(2n))}.
\]

Its2n equilateral triangles are(V_i,V_{i+1},U_i) and(V_i,U_{i-1},U_i). Every base link is(n,3,3,3).

In each case the U-ring is a regular unit-edge n-gon and supplies an explicit convex cap. Both completed surfaces have strictly positive defect2pi/n at every vertex. Their curvature is exactly4pi. Adjacent prescribed faces need not have been assumed noncoplanar during the exclusions.

All2n vertices of either belt lie on the sphere with center and squared radius

\[
c=(0,0,H/2),\qquad \rho^2=R^2+H^2/4,
\]

taking H=1 for the prism. This concludes the geometric and extendibility classification.

## 7. Verification and relation to earlier work

The exact enumeration and certificates are:

* `verify-general-local.py`: unrestricted local-link inequalities and exact hexagon-gap separation for n>=12;
* `verify_links11.py`: all16 unoriented/29 directed n=11 links;
* `verify_n11_gram.py`: both exact supporting-plane factor identities for the directed Q2 transitions;
* `../belt-study-12/verify_triangular_identities.py`: the symbolic triangular-fold and completion-plane identities;
* `../belt-study-12/certify_alternating45.py`: the n=12 square/pentagon compatibility exclusion;
* `certify_n11_q2_pair.py`: exact n=11 state comparisons, convex-root choice, and the gap between alpha34 and alpha35.

All listed checks passed. An [independent review](independent-review.md) audited the state ranges, transition uniqueness, parity, and completion-plane arguments and found no substantive gap.

Exploratory scripts are not used as certificates. The proof uses only necessary fragment conditions until the final explicit caps, and does not equate independent cone feasibility with admissibility.

For historical context, Grünbaum and Johnson's classification-free face-restriction paper treats the strict-face n>=12 and n=11 cases in Lemmas1–2 ([original paper](https://faculty.washington.edu/moishe/branko/BG53.Regular-faced%20polyh.pdf)). It was consulted as a structural cross-check, not used as a black-box completeness argument. The present proof supplies the coplanar boundary treatment, a direct fragment-level triangular obstruction, and an n=11 exclusion allowing unbounded future face sizes.
