# Changes from BI numeric 3.0 to BI certifier 4.0

- One shared contraction proof replaces the two duplicated implementations.
- Proof conditions are explicit exceptions, active even with Python optimization.
- Certificate rationals reject floats/booleans; triangle duplicates, label types and base labels are checked explicitly. Reversed intervals are rejected before rounding.
- Integer interval endpoints bypass general Fraction construction, using exact scaling by 2^bits.
- Exact integer dot products round outward once, giving tighter enclosures with fewer arithmetic objects.
- Exact chord constants are cached separately from precision-dependent intervals. The unused full chord-length table and dense Cartesian Jacobian are removed.
- Optional internal fold-angle computation, π enclosure and inverse-cosine bisection are removed. Existing nonflatness exclusions use chord/height witnesses.
- A single proposal module tries original-metric and BI-trace routes. Redundant internal certificate replays are removed; mandatory independent final replay is retained.
- Direct oriented face-adjacency propagation replaces the general graph matcher, removing NetworkX from runtime dependencies. Complete permutation/orbit comparison covers all 206 saved inputs.
- Dense Python and C++ control flow is formatted for review. Unused imports are removed. The numerical gradient algorithm and bundled geometry library are retained.
- Census 1.16 candidate-directory discovery and numerical pass ordering are supported.
- The standalone build can select the installed macOS Command Line Tools without modifying system settings; recent SDK libc++ headers are handled explicitly.
- Mathematical audit, adversarial arithmetic/certificate tests and corpus validation tools are included.

The two mathematical endpoint methods, exact intrinsic flip replay, complete patch-triangulation enumeration, strict convexity inequalities, original-face witnesses and conservative unresolved statuses are preserved. This release does not add exact-flatness acceptance to BI.
