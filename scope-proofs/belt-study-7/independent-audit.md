# Independent audit of the heptagon study

The local enumeration and all component-exclusion files were independently
read and their certificates executed. No substantive gap was identified.
The final assembly must cite these component arguments in dependency order.

## Completed checks

1. The angle enumeration is exhaustive without an empirical cutoff. The
   strict inequalities at the base give38 links (7,3,k),6 links (7,4,k),
   two links (7,5,k), and1+2+2 tetravalent links modulo reversal, totaling51.
   Expanding directions gives97. Independently ran the enumerator: PASS.
2. Read the generic Gram derivation and its sign conventions. The stated
   supported branch, signs before squaring, quadratic endpoints, and
   coplanar endpoints are consistent. Ran the state certificate: PASS.
3. The 333--Q1 bridge certificate considers both geometric roots for its
   necessary two-triangle completion. Its cutting-plane witnesses are
   actual vertices of the retained base or existing connector square.
   Ran the certificate: PASS, including q=0.
4. The all-pentagon trivalent component is excluded by two forced heptagons
   at a central pentagon vertex, whose total angle is71pi/35>2pi. The
   distinctness of those heptagons follows from proper face intersections.
5. The 334 study has exact finite strict-gap certificates and exact
   symbolic fixed-state identities. The forced pair of heptagons at a
   Q1-square corner has nonidentical outgoing rays and less than60 degrees
   of remaining face-angle budget. This is an obstruction to compatible
   neighboring completions, even though each forced face individually
   supports the original fragment. Its final odd-parity argument is
   explicitly conditional on the separately proved cross-component bridge
   exclusions.
6. Independently read `pentagon-proof.md` and `pentagon_certificate.py`.
   The P1--Q1 interval cover, exhaustive P1--Q2 root enclosure, and both
   necessary two-triangle completion roots are handled correctly. The
   pruning procedure need not assume uniqueness of the initial root.
   Independently ran the entire certificate: PASS (1813 rational boxes in
   the P1--Q2 portion).
7. At the fixed335 state, the P2--P2 completion test includes both geometric
   roots and both3+5 orders by reflection. The remaining shoulder is an
   actual vertex of the old edge pentagon and the new completing pentagon;
   its276-degree existing angle and gap strictly between60 and90 exclude
   further completion. The opposite root is rejected by a face plane
   cutting explicitly constructed existing vertices. This is a compatible
   completion test, not independent-cone reasoning.
8. The335 odd-state argument excludes trivalent K and hexagon blocks by
   exact state inequalities. Excluding both P1--P1 and P2--P2 then forces
   the period-three link pattern P1,P2,P2, incompatible with seven vertices.
   The334 odd-state and corner arguments similarly have all bridge
   dependencies supplied after the335 bridge exclusions are proved.

A minor mismatch between the initial333--Q1 gap bound stated in prose and
its weaker first certificate version was reported and corrected by the
author. The current script certifies the stated alpha11 bound and the
next-state estimate over the whole initial interval.

## Final assembly

Read the assembled `PROOF.md` after the component notes were complete.
The dependency order is acyclic: the335 bridge exclusions and finite
trivalent tests precede the333--334 bridge, and only then are the three
odd-cycle involutions applied separately. Intervening trivalent blocks do
not reconnect the separated families. In particular, the shared potential
hexagon bridge from333 toQ1 is excluded by the finite K6--Q1 test; Q-family
triangle--K interfaces retain only K7, outside the P2 interval and excluded
at P1. Pentagons and hexagons on square/pentagon blocks cannot change
their adjacent face-size label between endpoints, by reflection.

The displayed surviving coordinates, unit-edge heights, cosphere equations,
positive vertex defects and total curvature4pi are correct. No new gap or
overclaim was identified in the final assembly.

These arguments track base-edge states and use existing vertices as
support witnesses. They do not assume a simple outer ring or distinctness
of all outer vertices. At the tested vertices, additional remote
incidences must consume the same remaining positive angle budget.
