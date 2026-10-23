"""Exact metric construction, intrinsic flip replay, topology, and edge equations.

Extended from the local exclusion pilot. No catalog or special-case geometry.
"""

from fractions import Fraction
from itertools import combinations
from exact_intervals import I, sum_products
from field import chord_squared
from proof import require


def edge(a, b):
    return tuple(sorted((a, b)))


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def dot(a, b):
    return sum_products(zip(a, b))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def initial_metric(polygons):
    """Fan triangulations with exact unit-regular-polygon lengths.

    Chord constants are reused within this call only: interval precision may
    change between certificates. No unused all-pairs metric is constructed.
    """
    triangles, lengths, exact, constants = [], {}, {}, {}
    for face in polygons:
        size = len(face)
        for i in range(1, size - 1):
            triangle = [face[0], face[i], face[i + 1]]
            triangles.append(triangle)
            for a, b in combinations(triangle, 2):
                gap = abs(face.index(a) - face.index(b))
                key = (size, min(gap, size - gap))
                if key not in constants:
                    value = chord_squared(*key)
                    constants[key] = value, value.interval().sqrt()
                value, length = constants[key]
                e = edge(a, b)
                require(e not in exact or exact[e] == value, "Conflicting intrinsic edge lengths")
                exact[e], lengths[e] = value, length
    return triangles, lengths


def replay_flips(polygons, flips):
    faces, L = initial_metric(polygons)
    require(isinstance(flips, list), "Expected a flip list")
    require(
        all(isinstance(f, list) and len(f) == 4 and all(type(v) is int for v in f) for f in flips),
        "Invalid flip labels",
    )
    for a, b, c, d in flips:
        require(len({a, b, c, d}) == 4, "Certificate condition failed")
        adjacent = [i for i, f in enumerate(faces) if a in f and b in f]
        require(len(adjacent) == 2, "Certificate condition failed")
        require(
            {frozenset(faces[i]) for i in adjacent} == {frozenset((a, b, c)), frozenset((a, b, d))},
            "Certificate condition failed",
        )
        require(edge(c, d) not in L, "pilot does not support parallel intrinsic edges")
        length = L[edge(a, b)]
        xc = (L[edge(a, c)].square() + length.square() - L[edge(b, c)].square()) / (2 * length)
        xd = (L[edge(a, d)].square() + length.square() - L[edge(b, d)].square()) / (2 * length)
        yc = (L[edge(a, c)].square() - xc.square()).sqrt()
        yd = (L[edge(a, d)].square() - xd.square()).sqrt()
        require(yc.lo > 0 and yd.lo > 0, "Certificate condition failed")
        intersection = (yd * xc + yc * xd) / (yc + yd)
        require(
            intersection.lo > 0 and (length - intersection).lo > 0,
            "unfolded quad must be strictly convex",
        )
        L[edge(c, d)] = ((xc - xd).square() + (yc + yd).square()).sqrt()
        del L[edge(a, b)]
        faces[adjacent[0]] = [c, d, a]
        faces[adjacent[1]] = [d, c, b]
    return faces, L


def topology(faces, V):
    require(type(V) is int and V >= 4, "At least four vertices required")
    require(isinstance(faces, list) and bool(faces), "Expected triangles")
    require(
        all(
            isinstance(f, list)
            and len(f) == 3
            and len(set(f)) == 3
            and all(type(v) is int and 0 <= v < V for v in f)
            for f in faces
        ),
        "Invalid triangular face",
    )
    require(len({frozenset(f) for f in faces}) == len(faces), "Duplicate triangles")
    counts = {}
    for fi, f in enumerate(faces):
        require(len(set(f)) == 3, "Certificate condition failed")
        for a, b in combinations(f, 2):
            counts.setdefault(edge(a, b), []).append(fi)
    require(all(len(x) == 2 for x in counts.values()), "Certificate condition failed")
    require(V - len(counts) + len(faces) == 2, "Certificate condition failed")
    for v in range(V):
        link = {}
        for f in faces:
            if v in f:
                a, b = [x for x in f if x != v]
                link.setdefault(a, set()).add(b)
                link.setdefault(b, set()).add(a)
        require(link and all(len(x) == 2 for x in link.values()), "Certificate condition failed")
        reached = set()
        todo = [next(iter(link))]
        while todo:
            x = todo.pop()
            if x not in reached:
                reached.add(x)
                todo.extend(link[x] - reached)
        require(reached == set(link), "Certificate condition failed")
    # Connectivity of the whole complex.
    reached = {0}
    while True:
        new = reached | {v for f in faces if reached.intersection(f) for v in f}
        if new == reached:
            break
        reached = new
    require(reached == set(range(V)), "Certificate condition failed")


def setup(cert):
    polygons = cert["polygons"]
    validate_input(polygons)
    V = max(v for f in polygons for v in f) + 1
    original_triangles, _ = initial_metric(polygons)
    topology(original_triangles, V)
    original_edges = {}
    angles = [0] * V
    for f in polygons:
        for v in f:
            angles[v] += Fraction(len(f) - 2, len(f))
        for a, b in zip(f, f[1:] + f[:1]):
            original_edges[edge(a, b)] = original_edges.get(edge(a, b), 0) + 1
    require(
        all(n == 2 for n in original_edges.values()) and all(a < 2 for a in angles),
        "Certificate condition failed",
    )
    for f, g in combinations(polygons, 2):
        common = set(f) & set(g)
        require(len(common) <= 2, "Certificate condition failed")
        if len(common) == 2:
            for h in (f, g):
                require(
                    any(set((a, b)) == common for a, b in zip(h, h[1:] + h[:1])),
                    "Certificate condition failed",
                )
    faces, L = replay_flips(polygons, cert["flips"])
    topology(cert["triangles"], V)
    require(
        {tuple(sorted(f)) for f in faces} == {tuple(sorted(f)) for f in cert["triangles"]},
        "Certificate condition failed",
    )
    topology(faces, V)
    require(
        isinstance(cert["base"], list)
        and len(cert["base"]) == 3
        and all(type(v) is int and 0 <= v < V for v in cert["base"])
        and len(set(cert["base"])) == 3,
        "Invalid base triangle",
    )
    a, b, c = cert["base"]
    require(frozenset((a, b, c)) in {frozenset(f) for f in faces}, "Certificate condition failed")
    ell = L[edge(a, b)]
    x = (ell.square() + L[edge(a, c)].square() - L[edge(b, c)].square()) / (2 * ell)
    y = (L[edge(a, c)].square() - x.square()).sqrt()
    require(y.lo > 0, "Certificate condition failed")
    fixed = {a: [I(0), I(0), I(0)], b: [ell, I(0), I(0)], c: [x, y, I(0)]}
    unknown = [(v, k) for v in range(V) if v not in fixed for k in range(3)]
    equations = [e for e in sorted(L) if not set(e) <= set(fixed)]
    require(len(equations) == len(unknown) == 3 * V - 9, "Certificate condition failed")
    return V, faces, L, fixed, unknown, equations


def positions(V, fixed, unknown, values):
    p = [[I(0) for _ in range(3)] for _ in range(V)]
    for v, x in fixed.items():
        p[v] = x
    for (v, k), x in zip(unknown, values):
        p[v][k] = x
    return p


def validate_input(polygons):
    require(isinstance(polygons, list) and polygons, "Expected a nonempty face list")
    require(
        all(isinstance(f, list) and len(f) in (3, 4, 5, 6, 8, 10) for f in polygons),
        "Certificate condition failed",
    )
    require(
        all(type(v) is int and v >= 0 for f in polygons for v in f), "Certificate condition failed"
    )
    vertices = {v for f in polygons for v in f}
    require(
        vertices == set(range(len(vertices))), "Vertex labels must be contiguous, starting at zero"
    )
    require(all(len(f) == len(set(f)) for f in polygons), "Certificate condition failed")
    require(len({frozenset(f) for f in polygons}) == len(polygons), "Certificate condition failed")


def orient_polygons(polygons):
    """Choose consistent orientation without changing labels or polygon geometry."""
    adjacency = {}
    for i, f in enumerate(polygons):
        for a, b in zip(f, f[1:] + f[:1]):
            adjacency.setdefault(edge(a, b), []).append((i, 1 if a < b else -1))
    require(all(len(items) == 2 for items in adjacency.values()), "Certificate condition failed")
    neighbors = [[] for _ in polygons]
    for items in adjacency.values():
        (i, a), (j, b) = items
        neighbors[i].append((j, -a * b))
        neighbors[j].append((i, -a * b))
    sign = {0: 1}
    todo = [0]
    while todo:
        i = todo.pop()
        for j, relative in neighbors[i]:
            wanted = sign[i] * relative
            if j in sign:
                require(sign[j] == wanted, "Nonorientable surface")
            else:
                sign[j] = wanted
                todo.append(j)
    require(len(sign) == len(polygons), "Certificate condition failed")
    return [f[:] if sign[i] > 0 else f[:1] + f[:0:-1] for i, f in enumerate(polygons)]
