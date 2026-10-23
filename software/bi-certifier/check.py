#!/usr/bin/env python3
"""Independent proof replay. Standard Python only; no numerical solver, census, or catalog."""

import argparse, json
from itertools import combinations
from pathlib import Path
import exact_intervals as iv
import metric as m
from field import chord_squared

from proof import precision, root_box, contraction_image, require


def bounds(x):
    return x.exact_endpoints()


SCHEMA = "crp-bi-certifier/4"


def norm2(v):
    return sum((x.square() for x in v), iv.I(0))


def certify_metric_root(cert):
    """Banach/Krawczyk certificate for the exact, gauged metric equations.

    Constants are exact intrinsic lengths enclosed by metric.setup, not rounded
    replacement lengths. Uniform interval bounds cover their true values.
    """
    precision(cert["bits"])
    V, triangles, lengths, fixed, unknown, equations = m.setup(cert)
    n = len(unknown)
    center, X = root_box(cert, n)
    p = m.positions(V, fixed, unknown, X)
    pc = m.positions(V, fixed, unknown, center)
    Fc = [norm2(m.sub(pc[a], pc[b])) - lengths[a, b].square() for a, b in equations]
    index = {vk: i for i, vk in enumerate(unknown)}
    sparse = []
    for a, b in equations:
        diff = m.sub(p[a], p[b])
        sparse.append(
            {
                index[v, k]: 2 * sgn * diff[k]
                for v, sgn in ((a, 1), (b, -1))
                for k in range(3)
                if (v, k) in index
            }
        )
    image, contraction = contraction_image(cert, center, X, Fc, sparse)
    return (
        m.positions(V, fixed, unknown, image),
        triangles,
        dict(
            variables=n,
            equations=len(equations),
            contraction_bound=contraction,
            root_box_radius=cert["radius"],
        ),
    )


def convex_support(triangles, p):
    """Sufficient strict convexity test, with explicit unresolved weak signs.

    No tolerance converts a zero-containing interval to a weak inequality.
    Since every triangle is then an exposed hull facet, and these distinct
    facets form the checked closed connected sphere, they exhaust the hull
    boundary. Thus this strict branch requires no guessed face merging.
    """
    triangles = m.orient_polygons(triangles)
    volume6 = sum((m.dot(p[a], m.cross(p[b], p[c])) for a, b, c in triangles), iv.I(0))
    orientation = 1 if volume6.lo > 0 else -1 if volume6.hi < 0 else 0
    distinct = all(norm2(m.sub(p[a], p[b])).lo > 0 for a, b in combinations(range(len(p)), 2))
    unresolved = []
    strict_failures = []
    nondegenerate = True
    for i, (a, b, c) in enumerate(triangles):
        normal = m.cross(m.sub(p[b], p[a]), m.sub(p[c], p[a]))
        n2 = norm2(normal)
        if n2.lo <= 0:
            nondegenerate = False
            unresolved.append(dict(triangle=i, reason="Triangle nondegeneracy not certified"))
            continue
        if not orientation:
            continue
        for v in range(len(p)):
            if v in (a, b, c):
                continue
            value = orientation * m.dot(normal, m.sub(p[v], p[a]))
            if value.hi < 0:
                continue
            record = dict(
                triangle=i,
                vertices=[a, b, c],
                vertex=v,
                oriented_determinant=bounds(value),
                signed_plane_distance=bounds(value / n2.sqrt()),
            )
            if value.lo > 0:
                strict_failures.append(record)
            else:
                unresolved.append(record)
    certified = bool(
        orientation and distinct and nondegenerate and not unresolved and not strict_failures
    )
    return dict(
        certified=certified,
        method="oriented-strict-triangle-support",
        distinct_vertices_certified=distinct,
        nondegenerate_triangles_certified=nondegenerate,
        signed_six_volume=bounds(volume6),
        orientation=orientation,
        undecided_supports=unresolved,
        strictly_violated_supports=strict_failures,
        note=(
            "All triangle supports certified."
            if certified
            else "Convex identification remains unresolved. Coplanarity is allowed by the input "
            "contract, but zero-straddling support intervals do not prove weak support."
        ),
    )


def measure_original_faces(polygons, p):
    """Only disprove flat isometric regular faces; never certify flatness.

    A flat image of an original regular polygon must preserve every chord.
    Therefore a chord defect also proves a folded original disk even when
    its boundary vertices alone do not detect the fold.
    """
    records = []
    for fi, f in enumerate(polygons):
        a, b, c = f[:3]
        normal = m.cross(m.sub(p[b], p[a]), m.sub(p[c], p[a]))
        n2 = norm2(normal)
        heights = []
        if n2.lo > 0:
            heights = [
                iv.I(0) if v in (a, b, c) else m.dot(normal, m.sub(p[v], p[a])) / n2.sqrt()
                for v in f
            ]
        volume_defects = [
            dict(vertex=v, signed_plane_distance=bounds(h))
            for v, h in zip(f, heights)
            if not h.contains_zero()
        ]
        chord_defects = []
        for i, j in combinations(range(len(f)), 2):
            gap = norm2(m.sub(p[f[i]], p[f[j]])) - chord_squared(len(f), j - i).interval()
            if not gap.contains_zero():
                chord_defects.append(dict(vertices=[f[i], f[j]], squared_chord_defect=bounds(gap)))
        records.append(
            dict(
                face=fi,
                vertices=f,
                size=len(f),
                flatness_disproved_at_metric_root=bool(volume_defects or chord_defects),
                status=(
                    "NONFLAT_AT_METRIC_ROOT" if volume_defects or chord_defects else "NOT_DISPROVED"
                ),
                signed_distances_to_first_three_vertex_plane=[bounds(h) for h in heights],
                nonzero_height_witnesses=volume_defects,
                nonzero_chord_witnesses=chord_defects,
            )
        )
    return records


def classify(convex, faces):
    defect = any(f["flatness_disproved_at_metric_root"] for f in faces)
    if not convex:
        return "UNRESOLVED_CONVEXITY"
    return "EXCLUDED_NONFLAT" if defect else "OPEN_FOR_ALGEBRAIC"


def verify(cert):
    require(isinstance(cert, dict), "Expected a certificate object")
    if cert.get("method") == "hull-metric-interval":
        import hull_metric

        return hull_metric.verify(cert)
    require(
        cert.get("method") == "interval-root", "Only numerical interval certificates are supported"
    )
    p, triangles, root = certify_metric_root(cert)
    convexity = convex_support(triangles, p)
    faces = measure_original_faces(cert["polygons"], p)
    status = classify(convexity["certified"], faces)
    return dict(
        schema=SCHEMA,
        verified=True,
        status=status,
        root_certified=True,
        convex_metric_certified=convexity["certified"],
        excluded=status == "EXCLUDED_NONFLAT",
        exact_face_flatness_certified=False,
        crp_membership_certified=False,
        needs_algebraic_stage=status != "EXCLUDED_NONFLAT",
        reason={
            "EXCLUDED_NONFLAT": "Certified convex realization has a nonzero original-face defect.",
            "OPEN_FOR_ALGEBRAIC": "Convex realization certified; original-face flatness was not disproved. No CRP acceptance is asserted.",
            "UNRESOLVED_CONVEXITY": "Metric-equation root enclosed; convex identification not certified. No exclusion is authorized.",
        }[status],
        coordinate_frame=dict(
            base=cert["base"],
            description="base[0] at origin; base[1] on positive x; base[2] in xy plane with positive y",
        ),
        coordinate_intervals=[[bounds(x) for x in row] for row in p],
        enclosure_target=(
            "convex Alexandrov realization"
            if convexity["certified"]
            else "metric-equation root; convex identification unproved"
        ),
        root_proof=root,
        convexity=convexity,
        faces=faces,
        internal_folds=[],
        internal_fold_note="Exclusion uses original-face chord and height witnesses; no angle claims.",
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    report = verify(json.loads(args.certificate.read_text()))
    if args.output:
        with args.output.open("x") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in (
                    "status",
                    "root_certified",
                    "convex_metric_certified",
                    "excluded",
                    "reason",
                )
            },
            indent=2,
        )
    )
