# Convex-hull metric interval certificate

This endpoint certificate uses BI only to propose a numerical neighborhood. It does not certify the BI iteration history or recover algebraic coordinates. Its construction and proof are described here so that the mathematical argument, not merely successful checker execution, can be reviewed.

## The map being solved

Let T be the intrinsic triangulation reconstructed from the input regular polygons and the checked flip trace. Let its edge lengths be l_e. After a six-coordinate Euclidean gauge, the coordinate variable x has dimension n=3V-6, equal to the number of edges of T. Unlike the ordinary Cartesian root certificate, the lengths of the base triangle are not fixed coordinates: they are included among the equations.

The certificate supplies a rational center c and a box X. For each x in X, let P(x) be the convex hull of the labeled points p_v(x). We construct a continuous, finitely piecewise continuously differentiable function G on X whose components are squared intrinsic lengths of the marked edges of T on the boundary of P(x). The equation is

    F_e(x) = G_e(x) - l_e² = 0.

Thus a root is a convex realization by construction, rather than a spatial metric root whose convexity would need a subsequent decision about zero dihedrals.

## Checks defining the whole neighborhood

The checker first validates the input, the flips, T, and sphere topology as in the original interval checker. T is coarsened into polygonal patches, each a disk with all of its vertices on its boundary. The checker verifies that T restricts to a noncrossing triangulation of every patch, that patches form a polygonal sphere, that intersections are at most a common boundary edge, and that every labeled vertex belongs to the coarse graph and has degree at least three.

At the center, a proposed triangular convex hull is checked using strict outward support against every other vertex. It must triangulate exactly these patches. This supplies a known embedding of the coarse graph and its rotation system on a sphere.

Throughout the entire box X, directed interval inequalities then prove:

1. A fixed tetrahedral determinant is separated from zero, so the hull remains three-dimensional. The gauge remains nondegenerate with its prescribed signs.
2. Each vertex is strictly exposed. A fixed exposing vector is the sum of its incident outward center-triangle normals; its dot product with p_v(x)-p_w(x) is strictly positive for every w≠v.
3. Every coarse boundary edge ab stays an exposed hull edge. With third vertices c,d, the variable vector

       N(x) = (p_b-p_a) × (p_c-p_a) + (p_d-p_a) × (p_b-p_a)

   is identically perpendicular to p_b-p_a. One sign of N has strictly negative dot product with p_v-p_a for every v outside {a,b}. Hence its supporting plane meets the vertex set exactly at a,b. This uses a cross-product identity, not a numerical claim that a measured face determinant is zero.

The coarse graph is therefore embedded by the same exposed edges in the boundary of P(x) for every x in the connected box. All its vertices stay distinct and extreme. Its spherical rotation system cannot change continuously without edges or their directions colliding; such a collision would require a nonincident intersection of distinct hull edges, or three collinear hull vertices, both impossible for the certified extreme vertices. Consequently the complementary disks have the same labeled boundary cycles as at c. There are no other vertices in their interiors because every labeled hull vertex belongs to the coarse graph.

Each such hull patch may have further edges and flat diagonals. Whatever its actual facets, it admits a noncrossing triangulation using its boundary vertices. Therefore its actual triangulation is among the complete Catalan enumeration checked below. No particular diagonal is asserted to be flat or to be an actual hull edge.

## Interval derivatives across all diagonal choices

For each patch and each of its noncrossing triangulations, form the Euclidean triangle-edge lengths from x and develop its intrinsic disk into the plane. The dual graph of a polygon triangulation is a tree, so this uses only additions, products, division by positive squared edge lengths, and positive square roots. Interval automatic differentiation bounds both values and derivatives. Every square root used in differentiation is certified strictly positive throughout X.

For every development and every boundary edge, strict directed planar support tests against all other developed vertices prove that its boundary polygon is simple and strictly convex throughout X. Thus every marked T diagonal lies in this intrinsic patch and its length is the Euclidean length between the developed endpoints. Its squared length supplies the corresponding component of G. Coarse boundary edges use their squared spatial length directly.

For the actual hull triangulation these formulas give the actual intrinsic patch metric. If the hull changes triangulation, the change occurs through coplanar facets. The two developments then represent the same planar metric subdivision, and their marked lengths agree. Hence G is continuous across the cell boundaries. The candidate formulas are C¹ on a neighborhood of X by the strict domain checks. There are finitely many hull combinatorial pieces, and the interval Jacobian contains the derivatives of every candidate formula. The resulting continuous piecewise-C¹ map is Lipschitz, and its differences along any segment in X obey the corresponding interval mean-value bound. This does not require selecting an exact diagonal pattern at the root.

Different patches can change independently. It is sufficient to take the interval hull of the derivative rows over all triangulations of each patch. This overbounds every combination without enumerating the Cartesian product of the patches' triangulations. The implementation enumerates each patch completely; exceeding the resource cap returns unresolved rather than dropping choices.

## Root existence and coordinate enclosure

Let Y be the rational proposed inverse and let J(X) enclose all of the piece derivatives. The checker verifies

    ||I - Y J(X)||∞ < 1

and strict inclusion of

    K = c - Y F(c) + (I - Y J(X))(X-c)

in X. The value F(c) is evaluated using the separately verified strict center hull, with the exact input lengths enclosed outward.

The map H(x)=x-YF(x) is therefore a contraction mapping X into itself. Banach's fixed-point theorem gives its unique fixed point x*. The same bound at a differentiable piece implies that Y is invertible (I-YA has norm less than one for a square A), so F(x*)=0. The mean-value bound also places x* in K, which is the reported, usually much narrower, coordinate enclosure.

By construction P(x*) is convex. Its marked patch disks have the T edge lengths l_e. Each of their T triangles is therefore congruent to the corresponding input intrinsic triangle. Gluing these congruences across the checked marking gives an isometry from the input metric to the convex hull boundary. Alexandrov uniqueness identifies this realization with the intended convex Alexandrov realization up to Euclidean motion.

No step concludes that an original face or an internal diagonal is exactly flat. Coplanar original faces and changing artificial diagonals are allowed.

## Numerical exclusion and remaining scope

After this convex identification, an original-face vertex height or regular-polygon squared-chord defect whose interval excludes zero disproves a flat image of that original disk. If no such witness is found, the result is OPEN_FOR_ALGEBRAIC, never an acceptance.

The actual hull diagonal pattern is deliberately not selected at the root. This certificate therefore does not report angles at guessed diagonals. Both endpoint methods report height/chord witnesses; version 4 removes the optional fold-angle measurements.

This is a sufficient, independently replayable endpoint certificate, not proof-assistant formalization or a guarantee that the producer always finds a suitable box and patch decomposition. Intrinsic flips unsupported by the metric replay, non-disk patch proposals, failed strict neighborhood conditions, exhausted triangulation budgets, and timeouts remain unresolved. Fully flattened doubled-polygon Alexandrov limits are not covered by the full-dimensional hull certificate.

References for the underlying standard tools:

- Bobenko–Izmestiev, Alexandrov's theorem, weighted Delaunay triangulations, and mixed volumes, Theorem 1.1: https://www.numdam.org/item/10.5802/aif.2358.pdf
- The interval certificate uses the contraction-mapping argument explicitly above; its piecewise mean-value bound is justified by the finite smooth hull pieces and their continuity, rather than assuming a single differentiable Cartesian metric system across flips.
