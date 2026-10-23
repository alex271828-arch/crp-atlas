"""Interval certification of the piecewise-smooth convex-hull intrinsic metric.

See HULL_METRIC_PROOF.md. No coordinate recognition, symmetry, or zero tests.
"""

from itertools import combinations
from functools import lru_cache
import exact_intervals as iv
import metric as m
from proof import require, rational, root_box, contraction_image


class AD:
    __slots__ = ("v", "d")

    def __init__(self, value, derivative=None):
        self.v = iv.I(value)
        self.d = {} if derivative is None else derivative

    def __add__(self, other):
        b = asad(other)
        d = dict(self.d)
        for k, v in b.d.items():
            d[k] = d.get(k, iv.I(0)) + v
        return AD(self.v + b.v, d)

    __radd__ = __add__

    def __neg__(self):
        return AD(-self.v, {k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self + -asad(other)

    def __rsub__(self, other):
        return asad(other) + -self

    def __mul__(self, other):
        b = asad(other)
        d = {k: v * b.v for k, v in self.d.items()}
        for k, v in b.d.items():
            d[k] = d.get(k, iv.I(0)) + self.v * v
        return AD(self.v * b.v, d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = asad(other)
        value = self.v / b.v
        d = {k: v / b.v for k, v in self.d.items()}
        for k, v in b.d.items():
            d[k] = d.get(k, iv.I(0)) - value * v / b.v
        return AD(value, d)

    def square(self):
        return AD(self.v.square(), {k: 2 * self.v * v for k, v in self.d.items()})

    def sqrt(self):
        value = self.v.sqrt()
        require(value.lo > 0, "Differentiated square root is not bounded away from zero")
        return AD(value, {k: v / (2 * value) for k, v in self.d.items()})


def asad(value):
    return value if isinstance(value, AD) else AD(value)


def squared_distance(a, b):
    return sum(((x - y).square() for x, y in zip(a, b)), AD(0))


def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def det(p, a, b, c, d):
    return m.dot(m.cross(m.sub(p[b], p[a]), m.sub(p[c], p[a])), m.sub(p[d], p[a]))


@lru_cache(maxsize=256)
def triangulations(cycle):
    """All noncrossing triangulations of a cyclic polygon, exactly once."""
    if len(cycle) < 3:
        return ((),)
    result = []
    for k in range(1, len(cycle) - 1):
        for left in triangulations(cycle[: k + 1]):
            for right in triangulations(cycle[k:]):
                result.append(((cycle[0], cycle[k], cycle[-1]),) + left + right)
    return tuple(result)


def catalan(n):
    from math import comb

    return comb(2 * n, n) // (n + 1)


def gauge(base, V):
    a, b, c = base
    require(
        len(set(base)) == 3 and all(type(v) is int and 0 <= v < V for v in base),
        "Certificate condition failed",
    )
    zero = {(a, 0), (a, 1), (a, 2), (b, 1), (b, 2), (c, 2)}
    return [(v, k) for v in range(V) for k in range(3) if (v, k) not in zero]


def positions(V, unknown, values, ad=False):
    p = [[AD(0) if ad else iv.I(0) for _ in range(3)] for _ in range(V)]
    for i, ((v, k), x) in enumerate(zip(unknown, values)):
        p[v][k] = AD(x, {i: iv.I(1)}) if ad else x
    return p


def prepare(cert):
    V, triangles, lengths, _, _, _ = m.setup(cert)
    patches = cert["patches"]
    require(isinstance(patches, list) and patches, "Certificate condition failed")
    require(
        all(isinstance(f, list) and len(f) >= 3 and len(set(f)) == len(f) for f in patches),
        "Certificate condition failed",
    )
    require(
        all(type(v) is int and 0 <= v < V for f in patches for v in f),
        "Certificate condition failed",
    )
    require(len({frozenset(f) for f in patches}) == len(patches), "Certificate condition failed")
    require(set(v for f in patches for v in f) == set(range(V)), "A vertex was lost inside a patch")
    for f, g in combinations(patches, 2):
        common = set(f) & set(g)
        require(len(common) <= 2, "Patch boundaries intersect improperly")
        if len(common) == 2:
            for h in (f, g):
                require(
                    any(set((a, b)) == common for a, b in zip(h, h[1:] + h[:1])),
                    "Certificate condition failed",
                )
    patches = m.orient_polygons(patches)
    fan = [[f[0], f[i], f[i + 1]] for f in patches for i in range(1, len(f) - 1)]
    m.topology(fan, V)
    # Check exactly that the metric triangulation coarsens to these disks.
    owner = {}
    boundary_incidence = {}
    for i, f in enumerate(patches):
        require(catalan(len(f) - 2) <= 2000, "Patch triangulation budget exceeded")
        for a, b in zip(f, f[1:] + f[:1]):
            boundary_incidence.setdefault(m.edge(a, b), []).append(i)
    require(all(len(x) == 2 for x in boundary_incidence.values()), "Certificate condition failed")
    degrees = {v: 0 for v in range(V)}
    for a, b in boundary_incidence:
        degrees[a] += 1
        degrees[b] += 1
    require(all(d >= 3 for d in degrees.values()), "A coarse vertex has degree below three")
    for t in triangles:
        candidates = [i for i, f in enumerate(patches) if set(t) <= set(f)]
        require(len(candidates) == 1, "Metric triangle has no unique patch")
        owner[frozenset(t)] = candidates[0]
    possible = []
    for i, f in enumerate(patches):
        choices = triangulations(tuple(f))
        given = {frozenset(t) for t in triangles if owner[frozenset(t)] == i}
        require(
            any(given == {frozenset(t) for t in ts} for ts in choices),
            "Metric triangles do not triangulate patch",
        )
        possible.append(choices)
    equations = sorted(lengths)
    unknown = gauge(cert["base"], V)
    require(len(equations) == len(unknown) == 3 * V - 6, "Certificate condition failed")
    diagonal_rows = [[] for _ in patches]
    boundary_rows = []
    for row, e in enumerate(equations):
        if e in boundary_incidence:
            boundary_rows.append(row)
        else:
            choices = [i for i, f in enumerate(patches) if set(e) <= set(f)]
            require(len(choices) == 1, "Certificate condition failed")
            diagonal_rows[choices[0]].append(row)
    seed = cert["seed_triangles"]
    require(len({frozenset(t) for t in seed}) == len(seed), "Certificate condition failed")
    m.topology(seed, V)
    seed_choices = []
    for i, f in enumerate(patches):
        given = {frozenset(t) for t in seed if set(t) <= set(f)}
        selected = [ts for ts in possible[i] if given == {frozenset(t) for t in ts}]
        require(len(selected) == 1, "Seed hull does not refine the patch graph")
        seed_choices.append(selected[0])
    require(sum(len(ts) for ts in seed_choices) == len(seed), "Certificate condition failed")
    return dict(
        V=V,
        patches=patches,
        triangles=triangles,
        lengths=lengths,
        equations=equations,
        unknown=unknown,
        possible=possible,
        seed_choices=seed_choices,
        diagonal_rows=diagonal_rows,
        boundary_rows=boundary_rows,
        boundaries=boundary_incidence,
    )


def verify_geometry(cert, ctx, X):
    """Every point of X has a full-dimensional convex hull with fixed ridges."""
    import check

    V = ctx["V"]
    center = [iv.I(rational(x)) for x in cert["center"]]
    pc = positions(V, ctx["unknown"], center)
    p = positions(V, ctx["unknown"], X)
    seed = check.convex_support(cert["seed_triangles"], pc)
    require(seed["certified"], "Center seed hull is not rigorously strictly convex")
    a, b, c = cert["base"]
    require(p[b][0].lo > 0 and p[c][1].lo > 0, "Gauge orientation is not stable")
    tetra = cert["interior_tetrahedron"]
    require(
        len(tetra) == 4
        and len(set(tetra)) == 4
        and all(type(v) is int and 0 <= v < V for v in tetra),
        "Certificate condition failed",
    )
    require(
        not det(p, *tetra).contains_zero(),
        "Hull dimensionality is not certified throughout the box",
    )
    # The sum of incident outward seed normals strictly exposes each vertex.
    normals = [[iv.I(0) for _ in range(3)] for _ in range(V)]
    oriented_seed = m.orient_polygons(cert["seed_triangles"])
    for a, b, c in oriented_seed:
        normal = m.cross(m.sub(pc[b], pc[a]), m.sub(pc[c], pc[a]))
        normal = [seed["orientation"] * v for v in normal]
        for v in (a, b, c):
            normals[v] = m.add(normals[v], normal)
    for v in range(V):
        for w in range(V):
            if v != w:
                require(
                    m.dot(normals[v], m.sub(p[v], p[w])).lo > 0,
                    "Vertex exposure is not certified throughout the box",
                )
    # A variable normal, perpendicular to the ridge for every point of X,
    # supports exactly its two endpoints. Third vertices are just witnesses,
    # not asserted to form facets at every point of X.
    ridges = cert["ridge_witnesses"]
    require(len(ridges) == len(ctx["boundaries"]), "Certificate condition failed")
    seen = set()
    for a, b, c, d, sign in ridges:
        require(type(sign) is int and sign in (-1, 1), "Certificate condition failed")
        require(
            all(type(v) is int and 0 <= v < V for v in (a, b, c, d)), "Certificate condition failed"
        )
        require(len({a, b, c, d}) == 4, "Certificate condition failed")
        e = m.edge(a, b)
        require(e in ctx["boundaries"] and e not in seen, "Certificate condition failed")
        seen.add(e)
        i, j = ctx["boundaries"][e]
        require(
            (c in ctx["patches"][i] and d in ctx["patches"][j])
            or (d in ctx["patches"][i] and c in ctx["patches"][j]),
            "Certificate condition failed",
        )
        ab = m.sub(p[b], p[a])
        normal = m.add(m.cross(ab, m.sub(p[c], p[a])), m.cross(m.sub(p[d], p[a]), ab))
        normal = [sign * x for x in normal]
        for v in range(V):
            if v not in (a, b):
                require(
                    m.dot(normal, m.sub(p[v], p[a])).hi < 0,
                    "Coarse ridge exposure is not certified throughout the box",
                )
    return p, dict(
        seed_strict_convexity_certified=True,
        vertices_exposed_throughout_box=True,
        ridges_exposed_throughout_box=len(seen),
        full_dimension_throughout_box=True,
    )


def develop(triangles, edge_squared):
    """Develop a triangulated polygon disk using interval automatic derivatives."""
    a, b, c = triangles[0]
    ab2 = edge_squared[m.edge(a, b)]
    ell = ab2.sqrt()
    x = (ab2 + edge_squared[m.edge(a, c)] - edge_squared[m.edge(b, c)]) / (2 * ell)
    y = (edge_squared[m.edge(a, c)] - x.square()).sqrt()
    points = {a: [AD(0), AD(0)], b: [ell, AD(0)], c: [x, y]}
    pending = list(triangles[1:])
    while pending:
        progress = False
        for i, t in enumerate(pending):
            known = [v for v in t if v in points]
            if len(known) != 2:
                continue
            for k in range(3):
                a, b, c = t[k], t[(k + 1) % 3], t[(k + 2) % 3]
                if a in points and b in points and c not in points:
                    break
            else:
                raise AssertionError("Invalid triangle orientation")
            pa, pb = points[a], points[b]
            dx, dy = pb[0] - pa[0], pb[1] - pa[1]
            ab2 = edge_squared[m.edge(a, b)]
            tval = (ab2 + edge_squared[m.edge(a, c)] - edge_squared[m.edge(b, c)]) / (2 * ab2)
            hval = (edge_squared[m.edge(a, c)] / ab2 - tval.square()).sqrt()
            points[c] = [pa[0] + tval * dx - hval * dy, pa[1] + tval * dy + hval * dx]
            pending.pop(i)
            progress = True
            break
        require(progress, "Triangulated patch did not develop as a disk")
    return points


def convex_development(cycle, points):
    # Strict supporting oriented lines prove that the cyclic polygon is convex
    # and simple. The combinatorial noncrossing triangulations thus realize it.
    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
        direction = [points[b][k].v - points[a][k].v for k in range(2)]
        for v in cycle:
            if v not in (a, b):
                displacement = [points[v][k].v - points[a][k].v for k in range(2)]
                require(
                    cross2(direction, displacement).lo > 0,
                    "Intrinsic patch development is not uniformly strictly convex",
                )


def functions(ctx, values, all_pieces):
    """F(center) on its certified hull or interval Jacobian over all hull pieces."""
    p = positions(ctx["V"], ctx["unknown"], values, ad=True)
    edge_squared = {}

    def edge_value(e):
        if e not in edge_squared:
            edge_squared[e] = squared_distance(p[e[0]], p[e[1]])
        return edge_squared[e]

    n = len(ctx["equations"])
    rows = [None] * n
    for row in ctx["boundary_rows"]:
        rows[row] = edge_value(ctx["equations"][row])
    counts = []
    for i, cycle in enumerate(ctx["patches"]):
        # Triangular patches have no intrinsic diagonal rows.
        if not ctx["diagonal_rows"][i]:
            counts.append(1)
            continue
        for a, b in combinations(cycle, 2):
            edge_value(m.edge(a, b))
        choices = ctx["possible"][i] if all_pieces else [ctx["seed_choices"][i]]
        counts.append(len(choices))
        for triangulation in choices:
            points = develop(triangulation, edge_squared)
            if all_pieces:
                convex_development(cycle, points)
            for row in ctx["diagonal_rows"][i]:
                a, b = ctx["equations"][row]
                value = squared_distance(points[a], points[b])
                if rows[row] is None:
                    rows[row] = value
                else:
                    previous = rows[row]
                    keys = set(previous.d) | set(value.d)
                    derivative = {
                        k: iv.I(
                            min(previous.d.get(k, iv.I(0)).lo, value.d.get(k, iv.I(0)).lo),
                            max(previous.d.get(k, iv.I(0)).hi, value.d.get(k, iv.I(0)).hi),
                            raw=True,
                        )
                        for k in keys
                    }
                    rows[row] = AD(
                        iv.I(
                            min(previous.v.lo, value.v.lo), max(previous.v.hi, value.v.hi), raw=True
                        ),
                        derivative,
                    )
    require(all(row is not None for row in rows), "Certificate condition failed")
    F = [row.v - ctx["lengths"][e].square() for row, e in zip(rows, ctx["equations"])]
    J = [row.d for row in rows]
    return F, J, counts


def verify(cert):
    import check

    check.precision(cert["bits"])
    ctx = prepare(cert)
    n = len(ctx["unknown"])
    center, X = root_box(cert, n)
    _, geometry = verify_geometry(cert, ctx, X)
    Fc, _, _ = functions(ctx, center, False)
    _, J, counts = functions(ctx, X, True)
    image, contraction = contraction_image(cert, center, X, Fc, J)
    p = positions(ctx["V"], ctx["unknown"], image)
    faces = check.measure_original_faces(cert["polygons"], p)
    status = check.classify(True, faces)
    return dict(
        schema=check.SCHEMA,
        verified=True,
        status=status,
        method="hull-metric-interval",
        root_certified=True,
        convex_metric_certified=True,
        excluded=status == "EXCLUDED_NONFLAT",
        exact_face_flatness_certified=False,
        crp_membership_certified=False,
        needs_algebraic_stage=status != "EXCLUDED_NONFLAT",
        reason=(
            "Certified convex-hull metric realization has a nonzero original-face defect."
            if status == "EXCLUDED_NONFLAT"
            else "Convex-hull metric realization certified; flatness not disproved. Pass to algebraic certification."
        ),
        coordinate_frame=dict(
            base=cert["base"],
            description="base[0] at origin; base[1] on positive x; base[2] in xy plane with positive y; six gauge coordinates fixed",
        ),
        coordinate_intervals=[[check.bounds(x) for x in row] for row in p],
        enclosure_target="convex Alexandrov realization",
        root_proof=dict(
            method="piecewise-smooth convex-hull metric Krawczyk",
            variables=n,
            equations=n,
            contraction_bound=contraction,
            root_box_radius=cert["radius"],
            patch_sizes=list(map(len, ctx["patches"])),
            patch_triangulations_checked=counts,
        ),
        convexity=dict(
            certified=True,
            method="convex-hull metric with uniformly exposed coarse ridges",
            **geometry
        ),
        faces=faces,
        internal_folds=[],
        internal_fold_note="Actual facet diagonals are not selected by this enclosure; exclusion uses certified original-face height/chord defects, not guessed fold edges.",
    )
