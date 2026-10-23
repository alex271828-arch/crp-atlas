# Exhaustive link growth: invariant and constructive argument

## Scope

The target class consists of finite connected convex polyhedral surfaces with
regular unit-edge faces of sizes 3,4,5,6,8,10, meeting edge-to-edge, and strictly
positive curvature at every vertex. Coplanar neighboring faces are always permitted; zero-curvature vertices are not. Each pair of distinct faces
meets in at most one vertex or one common edge. The seed is a complete star of
the pass's designated vertex type. All vertex links belong to that pass's
remaining alphabet. Vertex labels of an existing fragment remain distinct.

These are the target assumptions, not claims about every generated candidate.
The argument below proves coverage of this class by the combinatorial search
provided its implementation follows the specified operations. It is a written
proof argument with executable tests, not a machine-checked formal proof and
not a final geometric certification of every candidate.

## Representation

A state is a set of consistently oriented, simple polygon face cycles. Its
face adjacency graph is connected: the seed has this property, and each added
face shares a boundary edge with the current complex. Each undirected edge
has at most two incident faces, which traverse it in opposite directions.

At vertex v a face (...,a,v,b,...) contributes the directed link arc a->b,
labelled by its face size. A partial link may be a disjoint union of paths.
A completed link is one cycle. Branching, several cycles, and a cycle plus
another component cannot embed injectively in one final cyclic link and are
rejected. A vertex with multiple paths is allowed, including when it appears
as an isolated contact between two face fans (a boundary pinch).

## The one-face extension lemma

Let P be a state contained injectively in a target Q, and let v have an
incomplete link. Choose any terminal node u of one of v's directed link paths.
The existing surface has a face using edge v->u. Q has exactly one other face
across this edge, traversing it u->v. This face is not already in P, since the
link path has no outgoing arc at u.

`next_faces` enumerates this missing face: its size is a letter in an allowed
completion of v's link. Each other corner is either an existing incomplete
vertex or a fresh vertex. An existing completed vertex cannot be incident to
a missing face of Q. All eligible old vertices and all fresh/reuse patterns
are enumerated, with fresh labels assigned in order of first occurrence.
Simple faces never require repeating a vertex within the same new face.

The edge, corner, and face-intersection checks hold for the Q face because it
belongs to Q. The resulting state remains a subcomplex of Q, even if it now
has disconnected fans at another vertex. The induced surface rotation may be tested, including at pinched vertices,
using the necessary sphere-embedding condition below.

## Whole-link completion lemma

Repeatedly apply the one-face lemma at the same v. Each step adds a missing
incident face. Positive curvature bounds a link to at most five faces for the
given face menu, so the process terminates with v's link completed. Thus every
choice of incomplete vertex has at least one enumerated complete-link child
that remains contained in Q. Pivot selection is an ordering heuristic, not a
coverage assumption. There is no alternate-pivot fallback.

The engine may pause between single-face additions. An inner agenda item
contains both P and the selected v. It becomes an outer item only when v's
link is complete. This is an operational subdivision of whole-link growth,
not a change to its enumeration semantics.

Mandatory-face closure may interleave the single-face steps. Each inserted
face is present in every target completion, so the state remains inside Q.
The selected pivot is preserved until its link is complete; closure may itself
complete it. No alternate-pivot fallback or conjectured repair is needed.

## Local-domain filter

Every existing path must occupy a contiguous directed run in one permitted
cyclic word. Distinct paths must occupy disjoint runs separated by at least
one missing face: without a missing face, their distinct endpoint edges would
have to be identified. The implementation enumerates all matching run
placements, both global orientations, and every cyclic position. Failure to
pack them is therefore a necessary obstruction. A cycle must itself equal
an allowed complete word.

The tests compare this packing implementation with a separate reference that
permutes entire paths, distributes explicit nonempty gaps, and matches the
result against rotated/reflected words.

## Other pruning rules

* Repeated vertices in a face, repeated directed edges, and improper pairwise
  face intersections violate the stated target class.
* All angle arithmetic is in units of three degrees. A completed word s has
  defect `120 - sum(ANGLE[n] for n in s)`. The sum of minimum feasible defects
  over all existing vertices may not exceed 240 (720 degrees). Every future
  vertex contributes additional positive defect. No floating tolerance enters
  this test.
* At most two path components fit a link of degree at most five: each
  component contains a face and needs a separating missing face. The cyclic
  ordering of two components is unique. Completing each local rotation thus
  gives the unique rotation system on the existing connected edge graph.
  A subgraph of a sphere must have V-E+R=2, where R counts rotation face walks.
  This test permits pinches; it does not require a manifold boundary.
* Each unfinished link has a joint row of missing words, one per gap. Region
  corners reference individual gaps. Cyclic propagation requires neighboring
  corners to agree on their shared missing polygon size. All gaps at a vertex
  retain their joint row; their curvature is counted only once. Removing rows
  without support is necessary, even though it may not establish global
  satisfiability. A missing polygon's proper boundary run has at most n-1
  edges. The existence-only run DP is equivalent to the former zero-cost DP.
* A complete missing face forced by those rows is added before branching
  and canonicalization. If it violates a necessary link or intersection
  condition, the original fragment is impossible. Every target completion
  contains each inserted face, so the original fragment and its forced
  closure have exactly the same completions. The selected inner pivot is
  preserved; forcing faces does not choose an alternate pivot.
* A closed state with connected cyclic links and consistent orientation has
  Euler characteristic 2, selecting a sphere.

## Conservative geometric pruning

`certified_geometry.py` uses fixed-point integer interval endpoints with
scale 2^52. Multiplication and division round down/up by integer floor and
ceiling. Square roots use integer square root and an explicit upper endpoint.
Face-angle cosines are rational or expressions in sqrt(2) and sqrt(5). No floating
geometric values, sampled lookup table or tolerance-based separation test
participates in rejection.

At a completed vertex, the convex spherical link has side lengths equal to
face angles. Degree-three links fix dihedral cosines. A degree-four link is
triangulated along a diagonal; a degree-five link along two diagonals. The
spherical cosine law and angle-addition identities give interval enclosures
for all transfers. Sines are nonnegative for angles in [0,pi]. A denominator
that may vanish returns an unknown cosine interval rather than excluding a
singular case. Convex subangle sums must not exceed pi. Inconsistent intervals
are the only geometric rejection mechanism.

Bounded caches retain local figure enclosures and proven exclusions, keyed
by the exact face word and integer interval endpoints. Quadrilateral keys
normalize cyclic rotation. Cached and uncached arithmetic have identical
outputs; cache reuse introduces no new geometric assumption. Each figure
cache is limited to 65,536 entries.

Local figure contraction subdivides one input cosine for degree-four links,
and two adjacent input cosines for degree-five links, into closed cells. Each
cell's output intervals must intersect the current bounds on every incident
edge. All surviving outputs are enclosed by their coordinatewise hull. Thus
local subdivision strengthens propagation before branching on the whole
fragment, without assuming independent realizability of those bounds.

When a spherical triangle has two identical algebraic side constants, its
base-angle cosine is evaluated as cot(a)*sqrt((1-cos(d))/(1+cos(d))). This
cancels the removable singularity at d=0. The diagonal satisfies
cos(d)>=2*cos(a)^2-1. For a=pi/2, the base-angle cosine is identically zero.
The identity is applied only to the same side constant, never to merely
overlapping intervals. The continuous extension encloses all nonsingular
triangles of a target's pointed convex vertex cone.

Angle sums and nonnegative differences use their monotone endpoint bounds in
cosine coordinates, rather than multiplying whole overlapping intervals in
an addition formula. Endpoints are still evaluated with outward-rounded
integer square roots and arithmetic. The sum is restricted to [0,pi], and
the difference to nonnegative angles; disjoint domains certify exclusion.

For each unfinished fan, `rim_geometry.py` develops the spherical endpoint
chord from its known internal dihedrals. The complementary chain of face
angles must span that same chord. The triangle inequality bounds its angular
length between max(0, largest side minus the other sides) and min(pi, sum of
sides). The union over supported completion words is necessary. Every fan at
a pinched vertex is checked separately, without assuming a single boundary
fan. Cosines of three-degree multiples use exact algebraic cos/sin(18-15)
and interval angle-addition recurrences, not floating trigonometry.

Each cyclic degree-four word also gets reusable individual edge-angle bounds.
The entire cosine interval [-1,1] is partitioned into 128 closed intervals.
A cell is discarded only if the same outward-rounded quadrilateral equations
contradict it; all remaining cells are retained by their hull. This is a cover
of a continuous domain, not sampling at 128 points. These bounds are intersected
at both endpoints of every shared edge before coupled propagation.

Regular polygons are developed from one fixed oriented face using enclosed
edge dihedral cosines. Convexity determines the fold sign. Coordinates of a
repeated vertex must intersect; each vertex must lie on the interior side of
every developed outward face plane. Unit vector components lie in [-1,1].
These give additional necessary interval tests, including for coplanar folds.

Free cosines start in [-1,1]. Proof-time budgets increase with the number of
completed links, to avoid expanding large unresolved parent fragments. Subdivision uses overlapping integer endpoint
boxes. Every unresolved box is retained, including on time or box-budget
exhaustion. Small boxes are shelved without claiming feasibility or exclusion,
while other boxes continue to be examined. When too many remain, their hull
is an outer enclosure. A parent cover is reused only under literal face-set
inclusion with the same vertex labels. Derived edge bounds can be inherited
because every parent's target realization satisfies them. An empty cover is
therefore a necessary impossibility certificate under these formulas.

Isomorphic states may retain different covers. Each cover encloses **all**
realizations of its fragment, so selecting either search representative is
safe. Covers are stored with inner and outer checkpoint items.

This is an implementation-level mathematical argument, not a separately
verified proof object or a machine-checked formal proof. Catalog and regression
tests exercise it but are not themselves an exhaustiveness proof. Final BI and algebraic certification are separate programs; Census journals
every closed candidate that survives these necessary filters.

## Isomorphism and the selected vertex

`canonical_code` selects starting face corners by isomorphism-invariant local
signatures. It tries all tied anchors in both orientations and traverses face
adjacencies deterministically, assigning vertex labels at first encounter.
The emitted list of labelled face cycles reconstructs the entire complex.
Equal codes therefore imply isomorphic complexes, and isomorphic complexes
admit corresponding anchors and produce the same minimum code.

Inner codes also include the canonical label of the selected vertex. Thus an
inner state is only merged with another state having the same rooted growth
problem. Outer states have no selected vertex. Their extensions are exhaustive
for any next pivot, so choosing a representative cannot lose a target.

Exact full codes are compared for deduplication; cryptographic hashes are
only human/machine identifiers in candidate journals, not search equality.

## Termination and induction

Let delta be the smallest allowed positive vertex defect, in three-degree units.
Gauss-Bonnet gives V <= floor(240/delta). The generator enforces this bound on
fresh labels, and all vertex links have degree at most five. Consequently the
number of faces and the number of labelled states are finite. Every face extension adds a face; the outer/inner bookkeeping performs only
finite additional work per extension. For any target containing the seed, the extension lemmas
preserve at least one path inside that target until it is closed. Each pass
therefore enumerates every target in its stated seed/alphabet class.

The standard seed needs no hidden old-vertex identifications: two distinct
faces in the star of a convex-polyhedral vertex already share that vertex;
additional intersections are restricted to their common incident edge.
Every target's designated minimum-type vertex supplies such a seed. Together
the passes cover the stipulated positive-curvature target class.

## What an exhausted run means

An empty agenda means that the combinatorial candidate enumeration exhausted.
It does not certify geometric realizability. Every surviving closed candidate
is written to `candidates_K.jsonl`; no numerical acceptance/rejection is run. `geometric_certification_complete` remains false. An exact final
geometric certifier is still needed for a fully certified classification.

Checkpoint writes are atomic. Inner states are saved with their pivot. All
sibling face extensions enter the agenda before an inner state is marked
expanded and checkpointed. A candidate is journalled before the checkpoint
can forget it. Replaying after a crash may revisit a journal entry but does
not lose it; the journal sink compares exact canonical codes.

## Preserved v1.15 alphabet and pass partition

The face menu is (3,4,5,6,8,10). Exact face angles in three-degree units are
20,30,36,40,45,48. A full turn is 120 units and total sphere curvature is
240 units. Strictly positive links still have degree at most five. Direct
enumeration of all cyclic words produces 63 rotation/reflection classes;
24 contain 8 or 10 and 39 are the original alphabet. Geometrically impossible
words are allowed in this combinatorial inventory and may be refuted later.

The smallest positive defect in the full alphabet is 2 units (6 degrees),
so the universal vertex bound is 120. Each pass uses its own remaining
alphabet to compute an equal or stronger bound. No vertex bound of 60 from
the old four-size alphabet is imposed on the new passes.

All new words precede the original 39 words, whose relative order is
unchanged. Every target containing an octagon or decagon contains at least
one new vertex word and hence belongs to the first block. Every target
without these faces belongs to the second block. Choosing the earliest
occurring vertex word partitions targets across the passes; seed choices
within a pass may still generate isomorphic candidates, which are deduplicated.

The old pass k is now k+24, with the same permitted alphabet. Checkpoints
are nevertheless incompatible because source fingerprints, the angle grid,
and the checkpoint version differ.

## v1.16 operational changes

The alphabet, pass order, exact code and interval formulas are preserved.
Catalogue lookup and numerical classification are removed from the runtime.
Boundary helpers for the superseded simple-boundary algorithm are removed;
the coupled rotation-system filter is the only boundary filter.

Mandatory-face closure has no arbitrary step cap. Each successful step adds
a distinct face, and its new local domains still satisfy the curvature bound.
Every locally admissible patch has at most 120 vertices and at most five
incident faces per vertex. Since faces have at least three corners, even a
noncompletable patch has at most 200 faces. A final spherical target has the
stronger bound of 182 faces from Euler's formula. Thus mandatory closure
terminates without assuming that the patch actually has a completion.
If a forced face is incompatible, the patch is rejected. Previously checked
old-old face intersections need not be checked again: only old-new pairs change.

An agenda item remains present until all its children and its candidate journal
entry have been computed successfully. Only then is it removed and its exact
code marked seen. An exception therefore cannot silently lose that item in
memory. Snapshots are written only between these transactions.

A journal record is newline-terminated, flushed and fsynced before any snapshot
can omit its generating item. A crash-torn last line is removed under the
pass lock and replayed from the older snapshot. A malformed complete line is
an error, not a record to skip. Atomic replacement prevents partially written
checkpoints from being read. These guarantees concern process interruption;
whole-filesystem loss or externally altered files are outside the model.
