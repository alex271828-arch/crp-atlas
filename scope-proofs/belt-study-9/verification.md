# Nonagon belt verification

The following checks were executed successfully with Python 3.10:

* links9-enumerate.py: exact enumeration of 22 marked links modulo reversal,
  40 directed links, with all size bounds derived from rational angle inequalities.
* links9-state-certificate.py: rational outward bounds for the state components,
  adjacent excluded face sizes, and a_10 < sqrt(2)R < a_11.
* quad-certificate.py: exact symbolic identities for the supported Q2 fold,
  both positive support signs, unit outgoing rays, and mixed-gap dot product
  -cos(2pi/9); exact inequalities for the Q1-Q1 exclusion.

An independent final audit of PROOF.md and quad-proof.md found no substantive
gap. It checked completeness of the component split, validity at n=9 of the
inherited Gram-sign assumptions and unique transition branches, coplanar
endpoints, absence of any maximum-face assumption, and the effect of possible
outer vertex identifications. Proper face intersections require the two
completing nonagons to be distinct in both the 388-degree and 370-degree
contradictions.

The explicit prism and antiprism coordinates in PROOF.md establish actual
extendibility and cosphericity. The computation verifies algebra and strict
comparisons; it is not a numerical search offered as proof of completeness.
