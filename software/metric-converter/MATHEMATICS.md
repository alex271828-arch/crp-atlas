# Exact coordinate conversion and verification

This is an executable certificate argument, not a proof-assistant formalization. The trusted mathematical operations are exact rational/algebraic arithmetic, positive-root isolation, and the geometric checks below. The producer and all saved source flags are untrusted.

## 1. Algebra and signs

The primitive polynomial is checked irreducible over Q, and its rational interval isolates exactly one real root (with a special exact rational-root case). The retained arithmetic from algebraic radii 3.0.11 normalizes coefficients modulo that polynomial. The stage-four override additionally certifies every embedding of a positive radical constant: SymPy may propose a field representation numerically, but the checker requires its square to equal the exact radicand power and its sign to be strictly positive. Nested roots are handled recursively. A different conjugate cannot be accepted merely because it satisfies the same minimal polynomial. This strengthens the inherited embedding routine; its original vendored source remains unchanged. Any additional regular-polygon constants are explicit real algebraic expressions; the coefficient field includes them and need not equal the radius field Q(alpha). Each parsed square root must have a proved nonnegative radicand, preventing non-real entries from entering a supposedly real coordinate system. Additional generators rᵢ satisfy rᵢ²=dᵢ>0, with the positive real square root selected. Sparse multiplication reduces repeated generators using these identities.

A coefficientwise zero is an exact equality. Sign claims require outward rational interval separation from zero; uncertainty is failure. The checker does not need to assume the supplied radicals are independent: redundant radicals may obstruct recognition of zero, but cannot turn a false equality into a proved one. Independent square-class recognition in the producer is only a compactness aid.

Arithmetic strings use a restricted AST parser, not eval or unrestricted parsing. Divisors are checked before cancellation. All proof conditions raise explicit failures and are active under Python optimization.

## 2. A radial coordinate proposal

Let Rᵢ be a certified candidate squared radius and dᵢⱼ² a prescribed original-face chord. Adjacent rays have scalar product gᵢⱼ=(Rᵢ+Rⱼ−dᵢⱼ²)/2. The first radial triangle gives a 3×3 Gram matrix G, whose leading principal minors must be positive.

For an oriented triangle (i,j,k), write u=qᵢ, v=qⱼ, a=Rᵢ, b=Rⱼ, c=gᵢⱼ, s=gᵢₖ, t=gⱼₖ and D=ab−c²>0. Set x=(bs−ct)/D and y=(at−cs)/D. In coordinates with metric G, n=G⁻¹(u×v) is orthogonal to u and v, and nᵀGn=D/det(G). The proposed new vertex is

qₖ = xu + yv + sqrt((Rₖ−xs−yt)det(G)/D) n.

The positive root realizes the chosen positive radial orientation. The first basis triangle fixes the overall handedness. This development is a proposal only: already known vertices are not accepted merely because a propagation path reached them. The independent global checker verifies the complete output, including every face chord and every sphere membership.

An existing Gram-coordinate representation bypasses this proposal construction but bypasses none of the global checks. Its center may be interior, on the boundary or outside; no unsupported interior-center assumption is applied to that route.

## 3. Independent geometric obligations

The checker does not import the reconstruction module. It verifies:

- Integer labels, the exact sphere partition, permissible regular face sizes, a connected closed combinatorial sphere, cyclic vertex links and strictly positive intrinsic curvature.
- Symmetry and positive definiteness of G by its leading principal minors.
- Distinctness of every pair of labeled vertices, by an exactly nonzero coordinate difference. Positive definiteness of G then implies positive squared separation.
- Every pairwise chord of every original regular polygon, with unit side length.
- Exact planarity of every face, and a nontrivial supporting halfspace containing every vertex. Nontrivial support also excludes a globally planar realization.
- Complete coverage of each coplanar supporting facet. Coplanar faces are oriented consistently; internal directed edges cancel; the remaining boundary must be a single strictly convex cycle through exactly all vertices on that supporting plane.
- Exact squared distance of every vertex from the recorded common center, equal to its assigned algebraic squared radius.

The facet argument follows the coordinate verifier in algebraic radii 3.0.11. Full regular-face distances imply each face is a simple convex regular polygon in its own plane. Summing their positively oriented boundaries cancels internal edges, leaving the strictly convex facet cycle. By planar winding, the polygons cover its interior once, without gaps or overlapping interiors. Positive intrinsic curvature rules out valid original vertices in a maximal facet's interior or in the relative interior of its boundary edges; the explicit boundary-vertex test enforces that condition here.

These complete supporting facets form the checked closed spherical surface and exhaust the convex hull boundary. The result is a full-dimensional convex realization with the prescribed regular flat faces. Coplanar neighboring original faces remain allowed.

## 4. Cartesian conversion

The producer factors G=L D Lᵀ with L unit lower triangular and D strictly positive diagonal, then sets B=sqrt(D)Lᵀ. The checker independently establishes BᵀB=G by exact sparse arithmetic and checks each exported Cartesian point equals Bqᵢ, including the center. Consequently every proved metric-coordinate property transfers to the ordinary Cartesian representation. This requires no approximate orthogonalization.

## 5. Representation claims

Printed formulas must be the deterministic rendering of the verified sparse coefficients and positive radicals. Coefficient factorizations are accepted only after exact equality checks. These checks prove equivalence of representations; they do not prove globally optimal simplicity. Neither field minimality nor expression minimality is asserted.

Every resource limit, unsupported algebraic situation, unresolved sign or failed proof obligation yields UNVERIFIED. Existing packages and their proof records are preserved.
