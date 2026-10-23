# Exact prism and antiprism families

For each integer n ≥ 3 let t be the largest real root of U_(n−1)(X), where U is the Chebyshev polynomial of the second kind. Its roots are cos(kπ/n), 1 ≤ k ≤ n−1; thus t = cos(π/n). Put s = √(1−t²)>0 and R=1/(2s). Use the polynomial recurrences displayed in the Atlas. They define algebraic coordinates, including a specific real embedding, without reliance on decimal evaluations.

For i=0,…,n−1 the lower vertices are
(R T_(2i)(t), R s U_(2i−1)(t), 0).
The upper vertices are
(R T_(2i+d)(t), R s U_(2i+d−1)(t), H).
For prisms d=0,H=1. For antiprisms d=1,H=√(1−1/[2(1+t)])>0.
The identity sin(kθ)=sin(θ) U_(k−1)(cos θ), with U_−1=0, identifies the two regular rings. Adjacent ring distances are 2R sin(π/n)=1. Prism sides are unit squares. Antiprism cross-edge squared distances are 2R²(1−t)+H²=1, so side triangles are equilateral.

## Supporting planes and convexity

Both ring planes support the convex hull. Prism side support planes are the supporting lines of the regular base polygon extended vertically; their intersections with the two rings give precisely the indicated four vertices.

For an antiprism set θ=π/n and take the lower edge with endpoint arguments 0 and 2θ. Let u=(cos θ,sin θ). Its endpoints have u-projection Rt, while the upper vertex at argument θ has projection R. The plane
u·(x,y) − R(1−t) z/H = Rt
contains this triangle. Every lower-ring projection is ≤ Rt, with equality precisely at that edge's endpoints; every upper-ring projection is ≤ R, with equality precisely at the indicated upper vertex. All vertices therefore lie in the same closed halfspace, and only the three stated vertices lie in the plane. Rotations and the symmetry interchanging upper and lower rings give all 2n triangle support planes. Together with the two bases these are exactly the hull faces, with the stated incidence. Thus these formulas produce convex regular-faced closed polyhedra for every n ≥ 3.

The common center is (0,0,H/2); all vertex squared radii are R²+H²/4. At each prism vertex the angle sum is (n−2)π/n+π=2π−2π/n. At each antiprism vertex it is (n−2)π/n+3π/3=2π−2π/n. Both have positive angular defect 2π/n.

The rendering limit n ≤ 256 is solely a browser display limit. The formulas and proof hold for all integers n ≥ 3. Some family members also occur in the finite census: the family views do not assert additional disjoint classes or a new count.

## Exhaustive scope

The finite census uses face sizes 3,4,5,6,8,10. The accompanying project reference proofs for sizes 7,9,and ≥11 give the full-belt restrictions used to reduce larger faces to these two families. Their existing proof and audit documents are copied into scope-proofs for reference. This Atlas is a presentation of the certified output and these family formulas; it does not rerun or replace the census completeness proof or its belt proofs.
