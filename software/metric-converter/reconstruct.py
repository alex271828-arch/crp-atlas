"""Untrusted exact proposal: glue positive radial cones, then verify globally."""

from collections import defaultdict, deque
import sympy as sp
from algebra import Algebra, require, chord, topology


def oriented(faces):
    faces = [list(f) for f in faces]
    edges = defaultdict(list)
    for i, f in enumerate(faces):
        for a, b in zip(f, f[1:] + f[:1]):
            edges[tuple(sorted((a, b)))].append((i, a, b))
    signs = {0: 1}
    todo = [0]
    while todo:
        i = todo.pop()
        f = faces[i]
        for a, b in zip(f, f[1:] + f[:1]):
            for j, c, d in edges[tuple(sorted((a, b)))]:
                if j == i:
                    continue
                s = signs[i] * (-1 if (a, b) == (c, d) else 1)
                if j in signs:
                    require(signs[j] == s, "Nonorientable face complex")
                else:
                    signs[j] = s
                    todo.append(j)
    require(len(signs) == len(faces), "Disconnected face complex")
    return [f if signs[i] == 1 else list(reversed(f)) for i, f in enumerate(faces)]


def reconstruct(case, saved):
    faces, V = topology(case, saved)
    A = Algebra(saved)
    F = A.field
    radii = [A.read_base(r) for r in saved["squared_radii"]]
    for r in radii:
        F.positive(r)
    rs = {v: r for group, r in zip(saved["groups"], radii) for v in group}
    if saved["verification_method"] == "exact_gram_geometry":
        geom = saved["geometric_certificate"]
        G = sp.Matrix([[A.read_base(x) for x in row] for row in geom["metric"]])
        points = [[A.base(A.read_base(x)) for x in row] for row in geom["coordinates"]]
        center = [A.base(A.read_base(x)) for x in geom["center"]]
        return A, G, points, center, {"route": "saved exact metric coordinates"}
    require(saved["verification_method"] == "radial", "Unsupported source method")
    faces = oriented(faces)
    triangles = []
    lengths = {}
    for face in faces:
        for i, u in enumerate(face):
            for j, v in enumerate(face):
                if i < j:
                    lengths[tuple(sorted((u, v)))] = F.normal(chord(len(face), j - i))
        triangles.extend((face[0], face[k], face[k + 1]) for k in range(1, len(face) - 1))

    def gram(i, j):
        return rs[i] if i == j else F.normal((rs[i] + rs[j] - lengths[tuple(sorted((i, j)))]) / 2)

    base = triangles[0]
    G = sp.Matrix([[gram(i, j) for j in base] for i in base])
    for k in range(1, 4):
        F.positive(F.normal(G[:k, :k].det()))
    determinant = F.normal(G.det())
    inverse = G.inv().applyfunc(F.normal)
    points = {v: [A.base(int(i == j)) for j in range(3)] for i, v in enumerate(base)}
    inv = [[A.base(x) for x in row] for row in inverse.tolist()]
    remaining = list(triangles)
    while remaining:
        changed = False
        pending = []
        for tri in remaining:
            missing = [v for v in tri if v not in points]
            if not missing:
                continue
            if len(missing) != 1:
                pending.append(tri)
                continue
            k = missing[0]
            at = tri.index(k)
            i, j = tri[(at + 1) % 3], tri[(at + 2) % 3]
            a, b, c = rs[i], rs[j], gram(i, j)
            u, v = gram(i, k), gram(j, k)
            D = F.normal(a * b - c * c)
            F.positive(D)
            x = F.normal((b * u - c * v) / D)
            y = F.normal((a * v - c * u) / D)
            h = F.normal((rs[k] - x * u - y * v) * determinant / D)
            F.positive(h)
            normal = A.matvec(inv, A.cross(points[i], points[j]))
            root = A.sqrt(h)
            points[k] = [
                A.total([A.scale(pi, x), A.scale(pj, y), A.mul(root, n)])
                for pi, pj, n in zip(points[i], points[j], normal)
            ]
            changed = True
        require(changed or not pending, "Could not propagate the radial triangulation")
        remaining = pending
    require(len(points) == V, "Incomplete coordinate reconstruction")
    return (
        A,
        G,
        [points[v] for v in range(V)],
        [{}, {}, {}],
        {"route": "positive radial cone development", "basis_vertices": list(base)},
    )


def cartesian(A, G, points, center):
    """G=L D L^T, so B=sqrt(D)L^T and B^T B=G."""
    F = A.field
    L = sp.eye(3)
    D = []
    for i in range(3):
        d = F.normal(G[i, i] - sum(L[i, k] ** 2 * D[k] for k in range(i)))
        F.positive(d)
        D.append(d)
        for j in range(i + 1, 3):
            L[j, i] = F.normal((G[j, i] - sum(L[j, k] * L[i, k] * D[k] for k in range(i))) / d)
    B = [[A.scale(A.sqrt(D[i]), L[j, i]) for j in range(3)] for i in range(3)]
    return B, [A.matvec(B, p) for p in points], A.matvec(B, center)
