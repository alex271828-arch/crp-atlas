# BI certifier 4.0: orbit constraints and radial gradient descent

The input is the same closed regular-polygon metric as in v2. No known coordinates or algebraic recognition are used.

## Symmetry and unknowns

`symmetry.py` enumerates the face-complex automorphisms by trying every image of a directed seed-face edge with either orientation, then propagating across uniquely paired face edges. Global vertex consistency and bijectivity are checked. Both orientation-preserving and orientation-reversing automorphisms are retained. The output records every permutation and the vertex orbits.

For a convex Alexandrov realization, these intrinsic isometries are realized by Euclidean isometries. A finite symmetry group has an invariant interior point, obtained by averaging the orbit of an interior point. Radii measured from that point agree on each vertex orbit. This justifies restricting the full-body search to orbit-invariant radii. If the group fixes a line or a plane, the restriction does not select a unique origin. Different vertex orbits may also have the same radius; the optimizer does not force or minimize those additional coincidences.

For an orbit O_a of size m_a, let x_a = r_i² for i in O_a. Let S be the square of the common initial radius. The solver uses scaled coordinates

    q_a = sqrt(m_a) x_a / S.

The sqrt(m_a) factor gives the Euclidean metric induced from the full vector of squared vertex radii, restricted to orbit-constant vectors. At every trial step, a single floating-point value r_a = sqrt(S q_a / sqrt(m_a)) is copied to the entire orbit. Symmetry is enforced before evaluating curvature, not imposed afterward by averaging a computed solution.

## Objective and gradient

For the radial tetrahedral complex, kappa_i = 2 pi minus the sum of dihedral angles around radial edge i. The fixed intrinsic Gaussian curvature at surface vertex i is a different quantity. The target is kappa_i = 0 for every radial edge.

    E(q) = (1/2) sum_i kappa_i(q)^2.

Writing J = d kappa / d r, the analytic gradient is

    (grad_q E)_a = S / sqrt(m_a) * sum_{i in O_a} (J^T kappa)_i / (2 r_i).

Geometry Central names the routine returning J `buildHessian`, because it is the Hessian of the BI functional. Here it is used only to evaluate the displayed gradient. No Newton equation is solved by the radial optimizer.

Every trial direction is exactly the negative gradient in q coordinates. Alternating Barzilai–Borwein scalar step lengths provide an initial step-size estimate; Armijo backtracking then requires both strict decrease of E and the sufficient-decrease inequality. A candidate that rounds to the current q is rejected as stalled. Thus a nonzero stationary point, inadmissible descent direction, or floating-point stagnation does not become a false success or an endless sequence of unchanged accepted steps.

## Geometry, edge flips, and restarts

Each accepted state must have positive radii and positive radial-tetrahedron Gram determinants, with a numerical margin relative to their scale. Intrinsic edges are flipped to remove boundary concavity. The flip queue is bounded. If a trial fails, every trial flip is discarded by rebuilding and replaying the last accepted history. Only accepted flips appear in the numerical handoff; the interval checker replays their metric independently.

By default the common initial radius is 1.5 times the largest initial triangle circumradius. The initializer increases it by 1.5 as needed to obtain a valid starting state. If descent fails, common-radius factors 2.0 and 1.2 are tried, with the same initialization check. At most three starts are used, stopping on the first success. If all fail, the initialized attempt with smallest final energy is retained as an uncertified estimate. Energy is monotone within each attempt, not across restarts. All restart outcomes are saved. `--radius-starts 1` requests a single start.

Default convergence means max_i |kappa_i| <= 1e-9 radians. An iteration or time limit is a resource bound, not a proof of nonexistence. Gradient descent on this nonlinear feasible domain is not guaranteed to reach zero curvature.

## Certification boundary

The symmetry group, numerical gradients, floating-point convexity checks, and convergence status are not part of the interval proof. `--numerical-only` stops after these operations and makes no certification claim.

The full pipeline passes the numerical coordinates to the Cartesian and hull proposal stages in `produce.py`. The producer first proposes a certificate using the original input triangulation with no flips, then tries the accepted BI trace if necessary. Using the original metric is a certificate proposal, not an assumption that the numerically recovered faces are flat: the verifier checks the complete intrinsic metric and convex realization for either route. This also avoids interval dependency growth from long flip histories and intermediate parallel edges that the old verifier does not support. That separate stage still uses Newton refinement to choose a proposed interval center. The exact-rational verifier independently checks the resulting certificate. It can certify a convex metric realization and sometimes disprove original-face flatness; it never infers exact face flatness merely from an interval containing zero. The endpoint proof obligations are retained and audited in `MATHEMATICAL_AUDIT.md`; the common contraction and exact arithmetic implementation are cleaned up in version 4.

The saved BI radii and apex describe the raw BI frame. Cartesian refinement uses a different gauge and does not overwrite those radii. Consequently, the floating-point radius vector and independently certified coordinate intervals must not be conflated.
