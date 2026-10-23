"""Untrusted numerical seeds for interval convex-hull metric certification."""

from fractions import Fraction as Q
from itertools import combinations
import numpy as np
from scipy.spatial import ConvexHull, QhullError
import exact_intervals as iv
import metric as m
import hull_metric as hm
import check


def patches_from_triangles(triangles, p):
    triangles = m.orient_polygons(triangles)
    normals = np.array([np.cross(p[b] - p[a], p[c] - p[a]) for a, b, c in triangles])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    parent = list(range(len(triangles)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    incidence = {}
    for i, t in enumerate(triangles):
        for a, b in combinations(t, 2):
            incidence.setdefault(m.edge(a, b), []).append(i)
    for i, j in incidence.values():
        if np.dot(normals[i], normals[j]) > 1 - 1e-10:
            parent[root(i)] = root(j)
    groups = {}
    for i, t in enumerate(triangles):
        groups.setdefault(root(i), []).append(t)
    patches = []
    for ts in groups.values():
        edge_counts = {}
        directions = {}
        vertices = set(v for t in ts for v in t)
        for a, b, c in ts:
            for v, w in ((a, b), (b, c), (c, a)):
                edge = m.edge(v, w)
                edge_counts[edge] = edge_counts.get(edge, 0) + 1
                directions[edge] = (v, w)
        boundary = [directions[e] for e, count in edge_counts.items() if count == 1]
        successor = dict(boundary)
        if len(successor) != len(boundary):
            raise ValueError("Proposed coarse patch is not a polygon disk")
        start = min(successor)
        cycle = [start]
        while successor[cycle[-1]] != start:
            v = successor[cycle[-1]]
            if v in cycle:
                raise ValueError("Proposed patch boundary repeats")
            cycle.append(v)
        if set(cycle) != vertices or len(cycle) != len(boundary):
            raise ValueError(
                "Proposed patch has an interior vertex or multiple boundary components"
            )
        patches.append(cycle)
    return patches


def propose(polygons, numerical, p, digits):
    patches = patches_from_triangles(numerical["triangles"], np.asarray(p))
    V = len(p)
    base = list(numerical["triangles"][0])
    unknown = hm.gauge(base, V)
    errors = []
    for exponent in (11, 10, 9):
        try:
            rng = np.random.default_rng(20260913)
            q = np.zeros((V, 3))
            for v, k in unknown:
                q[v, k] = p[v, k] + 10.0 ** (-exponent) * rng.uniform(-1, 1)
            seed = ConvexHull(q).simplices.tolist()
            center = [str(Q(float(q[v, k]))) for v, k in unknown]
            cert = dict(
                method="hull-metric-interval",
                bits=max(160, 4 * digits + 80),
                polygons=polygons,
                triangles=numerical["triangles"],
                flips=numerical["flips"],
                base=base,
                patches=patches,
                seed_triangles=seed,
                center=center,
                radius=str(Q(1, 10 ** (exponent - 3))),
            )
            check.precision(cert["bits"])
            ctx = hm.prepare(cert)
            # A tetrahedron robustly distinguishes the two spatial orientations.
            a, b, c = base
            normal = np.cross(q[b] - q[a], q[c] - q[a])
            d = int(np.argmax(np.abs((q - q[a]) @ normal)))
            cert["interior_tetrahedron"] = [a, b, c, d]
            ridge_witnesses = []
            for (a, b), (i, j) in ctx["boundaries"].items():
                best = None
                for c in ctx["patches"][i]:
                    for d in ctx["patches"][j]:
                        if len({a, b, c, d}) < 4:
                            continue
                        normal = np.cross(q[b] - q[a], q[c] - q[a]) + np.cross(
                            q[d] - q[a], q[b] - q[a]
                        )
                        values = np.array(
                            [np.dot(normal, q[v] - q[a]) for v in range(V) if v not in (a, b)]
                        )
                        for sign in (1, -1):
                            margin = -float(np.max(sign * values))
                            if best is None or margin > best[0]:
                                best = (margin, [a, b, c, d, sign])
                if best is None or best[0] <= 0:
                    raise ValueError("Numerical ridge support proposal failed")
                ridge_witnesses.append(best[1])
            cert["ridge_witnesses"] = ridge_witnesses
            C = [iv.I(Q(x)) for x in center]
            _, J, _ = hm.functions(ctx, C, False)
            dense = np.array(
                [
                    [
                        float(Q(row.get(j, iv.I(0)).lo + row.get(j, iv.I(0)).hi, 2 * iv.SCALE))
                        for j in range(len(unknown))
                    ]
                    for row in J
                ]
            )
            Y = np.linalg.inv(dense)
            cert["preconditioner"] = [[str(Q(float(x))) for x in row] for row in Y]
            return cert, check.verify(cert)
        except (
            AssertionError,
            ValueError,
            ArithmeticError,
            np.linalg.LinAlgError,
            QhullError,
        ) as error:
            errors.append(str(error))
    raise ValueError("Convex-hull metric proposal failed: " + "; ".join(errors))
