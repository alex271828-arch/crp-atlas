# Mathematical audit of Census v1.16

## Conclusion and scope

The retained filters are necessary conditions for the stated finite CRP°
class. The audit found no mathematical reason to remove or strengthen a
retained rejection rule. The new version preserves the spherical formulas
and the growth/pass semantics. Changes simplify implementation and improve
failure handling; they do not establish a new classification theorem.

The assurance remains an implementation-level proof argument backed by
independent arithmetic/combinatorial tests. This is not a proof-assistant
formalization or a separately replayable certificate for every discarded
search subtree. A successful regression against the old census cannot rule
out an error shared by both versions.

## Coverage obligations

1. **Finite target and pass partition.** Face angles in three-degree units
   are 20,30,36,40,45,48. Strict positive curvature gives degree at most five.
   All permitted cyclic words are enumerated modulo rotation/reflection.
   Their minimum defect is 2 and total sphere curvature is 240, so V≤120.
   Every target belongs to the pass of its earliest occurring word. No
   catalogue membership is used to choose seeds, menus or bounds.
2. **Injective seed.** Faces incident with the chosen vertex already share
   that vertex. The pairwise-intersection condition prevents an additional
   shared nonincident vertex or edge. Consecutive star faces share their
   single radial edge. The standard seed therefore needs no hidden vertex
   identifications.
3. **Extension.** At a selected unfinished radial edge the target has exactly
   one missing incident face. Its other vertices are fresh or existing open
   vertices; a closed vertex cannot receive another face. The generator tries
   all such choices, assigning fresh labels in first-occurrence order. It
   neither assumes a simple boundary nor changes pivot to escape a pinch.
4. **Local links.** A final cyclic link can contain disjoint directed paths
   only when they are separated by nonempty missing-face runs. Branching,
   multiple cycles or a cycle plus another component are impossible. The
   local support routine checks the entire joint packing, not each path alone.
5. **Partial rotations.** At most two nonempty path components fit a final
   degree-five link: each path requires at least one face and one gap.
   Their cyclic ordering is unique. Closing the missing adjacency in this
   rotation system fixes all complementary walks of the connected existing
   graph. A target sphere requires V−E+R=2. Repeated vertices on complementary
   walks are allowed; each corner refers to a particular gap, not merely
   to a vertex.
6. **Joint boundary words.** A completion supplies one row of missing words
   at each open vertex. Shared-edge face sizes must agree. Forward/backward
   propagation removes a row only if it has no supported closed cycle;
   filtering joint rows carries constraints between different gaps. Relaxing
   other dependencies can retain false positives but cannot reject a target.
7. **Runs, forcing and curvature.** Consecutive corners with one missing face
   belong to the same polygon. A proper run has fewer than n edges; an entire
   cycle can be that polygon only when its length is n. If the surviving
   possibilities force a complete n-vertex boundary, its closing face is
   mandatory. Every successful insertion strictly adds a face. Curvature is
   counted once per actual vertex, including pinches; future vertices only
   add positive curvature. Old-old intersections stay valid, so only pairs
   involving the newly forced face need rechecking.
8. **Deduplication.** The root-aware canonical code reconstructs the complete
   face complex. All tied isomorphism-invariant anchors and both orientations
   are tried. Building the edge map once per orientation leaves these choices
   unchanged. Inner codes include the pivot; outer pivots can be chosen freely
   because every pivot has an exhaustive completion procedure.

## Interval obligations

The stored interval [a,b] means [a/2^52,b/2^52]. Addition and negation are exact;
multiplication takes the extremal endpoint products and rounds outward;
division takes extremal endpoint quotients only when zero is absent from the
denominator; square handles a zero crossing explicitly; integer square root
rounds its upper endpoint upward. No binary floating-point geometric value
participates in a rejection. The added bounded sine cache keys on both exact
integer endpoints and returns the unchanged enclosure; it does not round
keys, discretize the domain, or identify merely overlapping intervals.

Regular-face constants are outward algebraic enclosures in √2 and √5.
An angle in [0,π] has the nonnegative sine √(1−cos²(angle)). For spherical
triangle side cosines x,y,z, the opposite angle cosine is
`(z−x*y)/(sqrt(1−x*x)*sqrt(1−y*y))`. Before any singular division the undivided
necessary equation is checked; a denominator interval containing zero returns
an unknown cosine bound. Identical fixed sides allow the algebraically
cancelled isosceles formula; overlapping interval values are not treated as
identical constants.

For sums, `cos(A+B)` decreases on 0≤A+B≤π. Feasibility is equivalent to
cos(A)+cos(B)≥0. For differences A−B≥0, feasibility is cos(A)≤cos(B).
The implementations use the monotone corner extrema, including −1 or +1 when
the constrained domain crosses the straight-angle or zero-angle boundary.
These bounds retain flat folds.

A spherical quadrilateral is split along one diagonal. A pentagon is split
into three triangles from two adjacent supplied angles. In a convex pointed
vertex cone the subangles have the prescribed nonnegative signs and add to
angles at most π. Coplanar adjacent original faces may give straight angles;
these are included. Distinct rays of such a target do not require arbitrary
signed square-root branches.

For incomplete links, the opposite chain's endpoint distance is at least
`max(0, largest side − sum(other sides))` and at most `min(π,sum(sides))` by
the spherical metric triangle inequality. Taking the union over completion
words is conservative, including when multiple partial fans are present.

Interval polygon development fixes the first face in a Euclidean frame and
uses the convex fold sign throughout. Repeated coordinates must intersect.
For every oriented face, all vertices must lie in its inner closed halfspace.
Clipping a unit-vector component to [−1,1] is a necessary condition, not a
normalization of a numerical approximation.

## Cover and resource invariants

- Each input box represents a continuous closed domain, not a sampled point.
- Bisections overlap at their midpoint, so no boundary value is lost.
- Local subdivision discards only contradicted cells and hulls all survivors.
- Time/box limits retain every queued or unresolved box. Narrow unresolved
  boxes remain unresolved rather than becoming feasibility assertions.
- Coalescing too many boxes takes an outer hull. If a variable is absent from
  one box, it is omitted from the hull restriction rather than invented.
- Parent boxes are inherited only under literal face-set inclusion with
  unchanged vertex labels. Every child realization restricts to a parent
  realization. Thus each retained cover contains all realizations of its
  fragment, and selecting one isomorphic representative remains safe.
- Only an empty remaining cover yields geometric exclusion.

## Changes deliberately avoided

No floating-point pruning, sampled angle table, new vertex bound, heuristic
face omission, symmetry restriction, midpoint branch selection or numerical
flatness test was introduced. Stronger algebraic pruning would increase the
proof burden and was outside this cleanup. Mirror-image enumeration remains
identified as before; an orientation-sensitive census would need a separate
canonicalization convention and validation.

## Remaining paper obligations

State the finite face-alphabet reduction and treatment of infinite families
explicitly. Decide whether to publish an implementation-level verified
computation or additionally build replayable subtree-exclusion certificates.
Keep enumeration exhaustion distinct from BI metric certification and exact
regular-face realizability. For the radial stage, distinguish its squared-radius
field from any claimed minimal coordinate realization field.
