#!/usr/bin/env python3
"""Independent coordinate certificate checker. Does not import the converter."""

import argparse
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import sympy as sp
from algebra import Algebra, require, chord, topology


def facet_patch(A, ids, faces, p, coplanar):
    """Positive polygon winding sums to the single convex facet boundary."""
    f = faces[ids[0]]
    normal = A.cross(A.vsub(p[f[1]], p[f[0]]), A.vsub(p[f[2]], p[f[0]]))
    occurrences = defaultdict(list)
    for i in ids:
        face = list(faces[i])
        a, b, c = [p[v] for v in face[:3]]
        orientation = A.sign(A.dot(normal, A.cross(A.vsub(b, a), A.vsub(c, a))))
        require(orientation != 0, "Degenerate coplanar polygon")
        if orientation < 0:
            face.reverse()
        for a, b in zip(face, face[1:] + face[:1]):
            occurrences[tuple(sorted((a, b)))].append((a, b))
    boundary = []
    for es in occurrences.values():
        if len(es) == 1:
            boundary.append(es[0])
        else:
            require(len(es) == 2 and es[0] == es[1][::-1], "Inconsistent facet edge cancellation")
    outgoing = {}
    incoming = {}
    for a, b in boundary:
        require(a not in outgoing and b not in incoming, "Branched facet boundary")
        outgoing[a] = b
        incoming[b] = a
    require(bool(outgoing) and set(outgoing) == set(incoming), "Open facet boundary")
    start = min(outgoing)
    v = start
    cycle = []
    while v not in cycle:
        cycle.append(v)
        v = outgoing[v]
    require(v == start and len(cycle) == len(boundary), "Multiple facet boundary components")
    require(
        len(cycle) >= 3 and set(cycle) == set(coplanar),
        "Facet boundary does not cover support vertices",
    )
    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
        for c in cycle:
            if c in (a, b):
                continue
            turn = A.dot(normal, A.cross(A.vsub(p[b], p[a]), A.vsub(p[c], p[a])))
            require(A.sign(turn) == 1, "Facet boundary is not strictly convex")
    return cycle


def verify(cert):
    require(cert.get("schema") == "crp-metric-coordinates/1", "Unsupported certificate schema")
    saved = cert["source_algebra"]
    case = cert["input"]
    faces, V = topology(case, saved)
    A = Algebra(saved)
    F = A.field
    for s in cert["positive_radicals"]:
        r = A.read_base(s)
        F.positive(r)
        A.rad.generators.append(r)
    require(len(A.rad.generators) <= 12, "Too many radical generators")

    def vector(v):
        require(isinstance(v, list) and len(v) == 3, "A vector needs three coordinates")
        return [A.decode(x) for x in v]

    def matrix(m):
        require(isinstance(m, list) and len(m) == 3, "A matrix needs three rows")
        return [vector(row) for row in m]

    G = sp.Matrix([[A.read_base(x) for x in row] for row in cert["metric"]])
    require(G.shape == (3, 3), "Metric must be 3 by 3")
    for i in range(3):
        for j in range(3):
            require(F.zero(G[i, j] - G[j, i]), "Metric is not symmetric")
    for k in range(1, 4):
        F.positive(F.normal(G[:k, :k].det()))
    gm = [[A.base(x) for x in row] for row in G.tolist()]
    p = [vector(v) for v in cert["metric_coordinates"]]
    center = vector(cert["metric_center"])
    require(len(p) == V, "Wrong vertex count")
    distances = {}
    for i, j in combinations(range(V), 2):
        delta = A.vsub(p[i], p[j])
        require(any(A.sign(x) != 0 for x in delta), "Coincident vertices or unresolved separation")
    facet_groups = defaultdict(list)
    regular_pairs = 0
    for fi, face in enumerate(faces):
        for i, j in combinations(range(len(face)), 2):
            pair = tuple(sorted((face[i], face[j])))
            if pair not in distances:
                distances[pair] = A.norm(A.vsub(p[pair[0]], p[pair[1]]), gm)
            d = distances[pair]
            want = A.base(chord(len(face), j - i))
            require(not A.sub(d, want), f"Incorrect regular-face chord on face {fi}")
            regular_pairs += 1
        a, b, c = [p[v] for v in face[:3]]
        normal = A.cross(A.vsub(b, a), A.vsub(c, a))
        sides = [A.sign(A.dot(normal, A.vsub(q, a))) for q in p]
        require(all(sides[v] == 0 for v in face), "Nonplanar face")
        nz = {s for s in sides if s}
        require(len(nz) == 1, "Face is not a supporting plane of a full-dimensional hull")
        coplanar = tuple(i for i, s in enumerate(sides) if s == 0)
        facet_groups[coplanar].append(fi)
    for vertices, ids in facet_groups.items():
        facet_patch(A, ids, faces, p, vertices)
    for group, r in zip(saved["groups"], saved["squared_radii"]):
        radius = A.read_base(r)
        F.positive(radius)
        for i in group:
            require(
                not A.sub(A.norm(A.vsub(p[i], center), gm), A.base(radius)),
                "Incorrect sphere membership",
            )
    B = matrix(cert["cartesian_transform"])
    cart = [vector(v) for v in cert["cartesian_coordinates"]]
    cc = vector(cert["cartesian_center"])
    require(len(cart) == V, "Wrong Cartesian vertex count")
    for i in range(3):
        for j in range(3):
            require(
                not A.sub(A.dot([r[i] for r in B], [r[j] for r in B]), gm[i][j]),
                "Cartesian transform is not an isometry",
            )
    for x, y in zip(p, cart):
        require(
            all(not z for z in A.vsub(A.matvec(B, x), y)),
            "Cartesian coordinates differ from metric coordinates",
        )
    require(all(not z for z in A.vsub(A.matvec(B, center), cc)), "Cartesian center mismatch")
    # Printed formulas are a deterministic rendering of verified sparse values.
    for name, values in [
        ("metric_coordinate_formulas", p),
        ("cartesian_coordinate_formulas", cart),
    ]:
        require(
            cert[name] == [[formula(A.encode(x)) for x in v] for v in values],
            "Printed formula mismatch",
        )
    return dict(
        status="CERTIFIED",
        verified=True,
        vertices=V,
        faces=len(faces),
        regular_face_pairs=regular_pairs,
        supporting_facets=len(facet_groups),
        exact_facet_tilings=True,
        exact_sphere_memberships=V,
        cartesian_isometry=True,
        printed_formulas_verified=True,
        positive_radicals=len(A.rad.generators),
        field_minimality_claimed=False,
        global_expression_optimality_claimed=False,
    )


def formula(element):
    if not element:
        return "0"
    terms = []
    for mask, c in sorted(element.items(), key=lambda x: int(x[0])):
        mask = int(mask)
        roots = [f"r{i}" for i in range(mask.bit_length()) if mask & (1 << i)]
        terms.append("(" + c + ")" + (" * " + " * ".join(roots) if roots else ""))
    return " + ".join(terms)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = verify(json.loads(args.certificate.read_text()))
    except Exception as exc:
        report = dict(
            status="UNVERIFIED", verified=False, reason=str(exc), exception=type(exc).__name__
        )
    if args.output:
        with args.output.open("x") as f:
            json.dump(report, f, indent=2)
            f.write("\n")
    print(json.dumps(report, indent=2))
    return 0 if report["verified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
