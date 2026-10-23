"""Exact local metric consequences of regular faces and body symmetries.

The equations are necessary conditions for a convex realization of the input.
No distance is estimated from numerical positions.  In particular, an inferred
square uses a combinatorial mirror and four proved unit distances, not a
four-cycle that merely looks square in a numerical drawing.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations

import sympy as sp


ALLOWED_FACE_SIZES = (3, 4, 5, 6, 8, 10)


def _edge(i, j):
    return tuple(sorted((int(i), int(j))))


def _zero(expression):
    return sp.simplify(expression) == 0


def _positive(expression):
    """An exact sign test; an undecided sign does not admit a metric patch."""
    expression = sp.simplify(expression)
    if expression.is_positive is not None:
        return bool(expression.is_positive)
    try:
        try:
            from .verify import ExactField, ProofFailure
        except ImportError:
            from verify import ExactField, ProofFailure
        return ExactField().sign(expression) == 1
    except (ValueError, TypeError, NotImplementedError, ProofFailure):
        return False


def squared_polygon_distance(n, step):
    """Exact squared chord of a unit-edge regular n-gon."""
    step = min(step % n, (-step) % n)
    if step == 0:
        return sp.S.Zero
    # Chebyshev avoids unrelated nested trigonometric radicals in chord ratios.
    cosine = sp.cos(sp.pi / n)
    return sp.simplify(sp.chebyshevu(step - 1, cosine) ** 2)


def _sphere_data(faces):
    """Check the hypotheses used when promoting automorphisms to mirrors."""
    edges = defaultdict(list)
    links = defaultdict(list)
    vertices = set()
    seen_faces = set()
    for face_id, face in enumerate(faces):
        if len(face) not in ALLOWED_FACE_SIZES or len(set(face)) != len(face):
            return None, "invalid regular face"
        if frozenset(face) in seen_faces:
            return None, "duplicate face vertex set"
        seen_faces.add(frozenset(face))
        for index, vertex in enumerate(face):
            vertices.add(vertex)
            following = face[(index + 1) % len(face)]
            edges[_edge(vertex, following)].append(face_id)
            links[vertex].append((face[index - 1], following))
    if len(vertices) < 4 or not faces or any(len(adjacent) != 2 for adjacent in edges.values()):
        return None, "input is not a closed edge-to-edge surface"
    if len(vertices) - len(edges) + len(faces) != 2:
        return None, "Euler characteristic is not two"
    face_graph = defaultdict(set)
    for left, right in edges.values():
        face_graph[left].add(right)
        face_graph[right].add(left)

    def connected(graph):
        if not graph:
            return False
        seen, pending = set(), [next(iter(graph))]
        while pending:
            current = pending.pop()
            if current not in seen:
                seen.add(current)
                pending.extend(graph[current])
        return len(seen) == len(graph)

    if len(face_graph) != len(faces) or not connected(face_graph):
        return None, "surface is disconnected"
    for vertex, pairs in links.items():
        link = defaultdict(list)
        for left, right in pairs:
            link[left].append(right)
            link[right].append(left)
        if any(len(neighbors) != 2 for neighbors in link.values()) or not connected(link):
            return None, "a vertex link is not a circle"
        angle_sum = sum(sp.Rational(len(face) - 2, len(face))
                        for face in faces if vertex in face)
        if angle_sum >= 2:
            return None, "a vertex lacks strictly positive intrinsic curvature"
    return dict(vertices=sorted(vertices), edge_faces=dict(edges)), None


def combinatorial_automorphisms(faces, *, max_automorphisms=256, max_flags=4096):
    """Enumerate cycle-preserving automorphisms by propagation from one face.

    The image and direction of one face determine the image across every edge
    of a connected closed surface.  This is a bounded flag enumeration, with no
    numerical geometry and no factorial permutation search.
    """
    faces = [tuple(map(int, face)) for face in faces]
    data, reason = _sphere_data(faces)
    if data is None:
        return {"permutations": [], "vertex_order": sorted(set(sum((list(f) for f in faces), []))),
                "complete": False, "reason": reason, "flags_examined": 0}
    vertices, edge_faces = data["vertices"], data["edge_faces"]
    size_counts = Counter(map(len, faces))
    initial = min(range(len(faces)), key=lambda i: (size_counts[len(faces[i])] * len(faces[i]), i))
    seed = faces[initial]
    found, seen = [], set()
    flags = 0

    def extend(target_id, shift, direction):
        vertex_map, inverse_vertices = {}, {}
        face_map, inverse_faces = {}, {}
        queue = deque()

        def assign_face(source_id, image_id, offset, sign):
            if source_id in face_map:
                return face_map[source_id] == image_id
            if image_id in inverse_faces or len(faces[source_id]) != len(faces[image_id]):
                return False
            source, target = faces[source_id], faces[image_id]
            for index, vertex in enumerate(source):
                image = target[(offset + sign * index) % len(target)]
                if ((vertex in vertex_map and vertex_map[vertex] != image)
                        or (image in inverse_vertices and inverse_vertices[image] != vertex)):
                    return False
                vertex_map[vertex] = image
                inverse_vertices[image] = vertex
            face_map[source_id] = image_id
            inverse_faces[image_id] = source_id
            queue.append(source_id)
            return True

        if not assign_face(initial, target_id, shift, direction):
            return None
        while queue:
            source_id = queue.popleft()
            image_id = face_map[source_id]
            source = faces[source_id]
            for index, left in enumerate(source):
                right = source[(index + 1) % len(source)]
                source_adjacent = edge_faces[_edge(left, right)]
                neighbor = source_adjacent[0] if source_adjacent[1] == source_id else source_adjacent[1]
                image_edge = _edge(vertex_map[left], vertex_map[right])
                if image_edge not in edge_faces:
                    return None
                image_adjacent = edge_faces[image_edge]
                if image_id not in image_adjacent:
                    return None
                image_neighbor = image_adjacent[0] if image_adjacent[1] == image_id else image_adjacent[1]
                if neighbor in face_map:
                    if face_map[neighbor] != image_neighbor:
                        return None
                    continue
                source_neighbor, target_neighbor = faces[neighbor], faces[image_neighbor]
                if len(source_neighbor) != len(target_neighbor):
                    return None
                sl, sr = source_neighbor.index(left), source_neighbor.index(right)
                tl, tr = target_neighbor.index(vertex_map[left]), target_neighbor.index(vertex_map[right])
                n = len(source_neighbor)
                source_sign = 1 if (sr - sl) % n == 1 else -1
                target_sign = 1 if (tr - tl) % n == 1 else -1
                sign = source_sign * target_sign
                offset = (tl - sign * sl) % n
                if not assign_face(neighbor, image_neighbor, offset, sign):
                    return None
        if len(face_map) != len(faces) or len(vertex_map) != len(vertices):
            return None
        return tuple(vertex_map[vertex] for vertex in vertices)

    for target_id, target in enumerate(faces):
        if len(target) != len(seed):
            continue
        for direction in (1, -1):
            for shift in range(len(seed)):
                if flags >= max_flags or len(found) >= max_automorphisms:
                    return {"permutations": found, "vertex_order": vertices,
                            "complete": False, "reason": "enumeration bound reached", "flags_examined": flags}
                flags += 1
                permutation = extend(target_id, shift, direction)
                if permutation is not None and permutation not in seen:
                    seen.add(permutation)
                    found.append(list(permutation))
    return {"permutations": found, "vertex_order": vertices,
            "complete": True, "reason": None, "flags_examined": flags}


def _canonical_cycle(cycle):
    cycle = tuple(cycle)
    candidates = []
    for oriented in (cycle, cycle[::-1]):
        candidates.extend(oriented[k:] + oriented[:k] for k in range(len(cycle)))
    return min(candidates)


def infer_reflection_squares(known_distances, reflection, *, vertices=None):
    """Force squares from a proved mirror; never infer one from a bare cycle.

    Two disjoint unit segments interchanged endpoint-for-endpoint by the mirror
    have endpoints at signed normal distance +/-1/2.  A unit connector cannot
    join opposite sides: that would make the segment midpoints coincide and
    identify two distinct vertices.  The connector is therefore on the same
    side, its midpoint distance is one, and the two remaining diagonals have
    squared length two.
    """
    pairs = reflection["swapped_pairs"]
    squares, seen = [], set()
    for (a, aa), (b, bb) in combinations(pairs, 2):
        if not (_zero(known_distances.get(_edge(a, aa), -1) - 1)
                and _zero(known_distances.get(_edge(b, bb), -1) - 1)):
            continue
        for left, right in ((b, bb), (bb, b)):
            if not (_zero(known_distances.get(_edge(a, left), -1) - 1)
                    and _zero(known_distances.get(_edge(aa, right), -1) - 1)):
                continue
            cycle = _canonical_cycle((a, left, right, aa))
            if cycle in seen:
                continue
            seen.add(cycle)
            diagonals = [_edge(a, right), _edge(aa, left)]
            # A contradicting already-known diagonal means this metric cannot
            # be such a nondegenerate square; do not overwrite it silently.
            if any(pair in known_distances and not _zero(known_distances[pair] - 2)
                   for pair in diagonals):
                continue
            squares.append({"vertices": list(cycle), "diagonals": diagonals,
                            "squared_diagonal": sp.Integer(2),
                            "unit_transposed_pairs": [[a, aa], [b, bb]],
                            "unit_connectors": [[a, left], [aa, right]],
                            "reflection_id": reflection["id"],
                            "proof": "unit mirror-transposed segments and unit same-side connectors force a nondegenerate square"})
    return squares


def _normalize_equation(expression, radius_symbols):
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    if numerator == 0:
        return sp.S.Zero
    polynomial = sp.Poly(sp.expand(numerator), *radius_symbols, extension=True)
    return sp.expand(polynomial.monic().as_expr())


def metric_patch_equation(anchors, known_distances, radius_at):
    """Return the exact four-point radius compatibility equation and metric.

    With q_i=p_i-p_0, G_ij=q_i.q_j and
    w_i=(s_0+|q_i|^2-s_i)/2=(-p_0).q_i, compatibility is
    w^T adj(G) w - s_0 det(G)=0.  All entries of G are already-known distances.
    """
    anchors = tuple(anchors)
    if len(anchors) != 4 or len(set(anchors)) != 4:
        return None
    if any(_edge(i, j) not in known_distances for i, j in combinations(anchors, 2)):
        return None
    base, *others = anchors
    radial = [sp.sympify(radius_at[v]) for v in anchors]
    lengths = [sp.sympify(known_distances[_edge(base, v)]) for v in others]
    gram = sp.Matrix(3, 3, lambda i, j: lengths[i] if i == j else
                     (lengths[i] + lengths[j] - known_distances[_edge(others[i], others[j])]) / 2)
    gram = gram.applyfunc(sp.simplify)
    minors = [sp.simplify(gram[:size, :size].det()) for size in (1, 2, 3)]
    if not all(_positive(value) for value in minors):
        return None
    w = sp.Matrix([(radial[0] + distance - radius) / 2
                   for distance, radius in zip(lengths, radial[1:])])
    adjugate = gram.adjugate()
    equation = sp.expand((w.T * adjugate * w)[0] - radial[0] * minors[-1])
    return {"anchors": list(anchors), "gram": gram.tolist(), "determinant": minors[-1],
            "positive_principal_minors": minors, "w": list(w), "equation": equation,
            "_adjugate": adjugate, "_gram_matrix": gram}


def enrich_system(model, positions=None, *, max_automorphisms=256,
                  max_flags=4096, max_cliques=4096, max_metric_patches=128):
    """Attach exact patch equations and their provenance to an encoder model.

    ``positions`` is accepted for caller compatibility but is never inspected.
    Truncation weakens the set of necessary equations and cannot manufacture a
    realization.  The original faces, equations, and sphere groups are retained.
    """
    faces = [list(map(int, face)) for face in model["faces"]]
    groups = [list(map(int, group)) for group in model["groups"]]
    radius_symbols = tuple(sp.sympify(x) for x in model.get("radii_symbols", model.get("radius_variables", [])))
    if len(radius_symbols) != len(groups):
        raise ValueError("Patch construction needs one squared-radius symbol per sphere group")
    radius_at = {vertex: radius for group, radius in zip(groups, radius_symbols) for vertex in group}
    vertices = sorted(set(sum(faces, [])))
    if (any(not group for group in groups) or sum(map(len, groups)) != len(vertices)
            or set(radius_at) != set(vertices)):
        raise ValueError("Sphere groups do not partition the face vertices")
    known, distance_sources, conflicts = {}, defaultdict(list), []

    def add_distance(pair, value, source):
        pair, value = _edge(*pair), sp.simplify(value)
        if pair in known and not _zero(known[pair] - value):
            conflicts.append({"vertices": list(pair), "first": known[pair], "second": value, "source": source})
            return False
        known[pair] = value
        distance_sources[pair].append(source)
        return True

    pair_faces = defaultdict(set)
    for face_id, face in enumerate(faces):
        for i, j in combinations(range(len(face)), 2):
            pair = _edge(face[i], face[j])
            pair_faces[pair].add(face_id)
            value = squared_polygon_distance(len(face), j - i)
            add_distance(pair, value, {"kind": "regular_face", "face": face_id,
                                      "face_size": len(face), "cyclic_step": min(j - i, len(face) - j + i)})

    symmetry = combinatorial_automorphisms(faces, max_automorphisms=max_automorphisms, max_flags=max_flags)
    reflections = []
    for permutation in symmetry["permutations"]:
        mapping = dict(zip(symmetry["vertex_order"], permutation))
        fixed = [vertex for vertex in vertices if mapping[vertex] == vertex]
        if len(fixed) < 3 or len(fixed) == len(vertices):
            continue
        if any(mapping[mapping[vertex]] != vertex for vertex in vertices):
            continue
        reflections.append({"id": len(reflections), "permutation": permutation,
                            "vertex_order": symmetry["vertex_order"], "mapping": mapping,
                            "fixed_vertices": fixed,
                            "swapped_pairs": [[v, mapping[v]] for v in vertices if v < mapping[v]],
                            "conditional_on": "existence of a closed convex regular-face realization with positive intrinsic curvature",
                            "proof": "the metric automorphism extends to an ambient isometry by convex metric uniqueness; three distinct extreme fixed vertices are noncollinear, so a nonidentity involution is a plane reflection"})

    squares, square_keys = [], set()
    for face_id, face in enumerate(faces):
        if len(face) == 4:
            square_keys.add(_canonical_cycle(face))
            squares.append({"vertices": list(_canonical_cycle(face)), "kind": "input_square", "face": face_id})
    inferred_squares = []
    for reflection in reflections:
        for square in infer_reflection_squares(known, reflection):
            key = tuple(square["vertices"])
            for pair in square["diagonals"]:
                add_distance(pair, sp.Integer(2), {"kind": "reflection_square", **square})
            if key not in square_keys:
                square_keys.add(key)
                square["kind"] = "reflection_square"
                squares.append(square)
                inferred_squares.append(square)

    pyramids = []
    for square in squares:
        boundary = square["vertices"]
        for apex in vertices:
            if apex not in boundary and all(_edge(apex, v) in known and _zero(known[_edge(apex, v)] - 1)
                                           for v in boundary):
                pyramids.append({"kind": "unit_square_pyramid", "apex": apex,
                                 "boundary": boundary, "square_source": square,
                                 "proof": "four unit apex distances to a proved unit square determine a nonflat square-pyramid metric"})

    neighbors = defaultdict(set)
    for left, right in known:
        neighbors[left].add(right)
        neighbors[right].add(left)

    def cliques():
        yielded = set()
        # The explicit template guarantees that a small pyramid patch is not
        # buried behind many coplanar cliques inside large regular faces.
        for pyramid in pyramids:
            anchor = (pyramid["apex"], *pyramid["boundary"][:3])
            key = tuple(sorted(anchor))
            if key not in yielded:
                yielded.add(key)
                yield anchor
        for a in vertices:
            for b in sorted(v for v in neighbors[a] if v > a):
                common = neighbors[a] & neighbors[b]
                for c in sorted(v for v in common if v > b):
                    for d in sorted(v for v in common & neighbors[c] if v > c):
                        key = (a, b, c, d)
                        if key not in yielded:
                            yielded.add(key)
                            yield key

    patch_equations, patch_sources, metric_patches = [], [], []
    equation_index = {}

    def add_equation(expression, source):
        equation = _normalize_equation(expression, radius_symbols)
        if equation == 0:
            return None
        if equation not in equation_index:
            equation_index[equation] = len(patch_equations)
            patch_equations.append(equation)
        index = equation_index[equation]
        patch_sources.append({"equation_index": index, "equation": equation, **source})
        return index

    examined, truncated = 0, False
    for anchors in cliques():
        if examined >= max_cliques or len(metric_patches) >= max_metric_patches:
            truncated = True
            break
        examined += 1
        common_faces = set.intersection(*(pair_faces[_edge(i, j)] for i, j in combinations(anchors, 2)))
        if common_faces:
            continue
        patch = metric_patch_equation(anchors, known, radius_at)
        if patch is None:
            continue
        patch_id = len(metric_patches)
        adjugate, gram = patch.pop("_adjugate"), patch.pop("_gram_matrix")
        patch["id"] = patch_id
        patch["distances"] = [{"vertices": [i, j], "squared_distance": known[_edge(i, j)],
                               "sources": distance_sources[_edge(i, j)]}
                              for i, j in combinations(anchors, 2)]
        patch["equation_index"] = add_equation(patch["equation"],
                                                {"kind": "four_point_metric_compatibility", "patch": patch_id,
                                                 "anchors": list(anchors),
                                                 "identity": "w^T adj(G) w - s_base det(G) = 0"})
        base, *others = anchors
        additional = []
        for vertex in sorted(set.intersection(*(neighbors[v] for v in anchors)) - set(anchors)):
            length = known[_edge(base, vertex)]
            h = sp.Matrix([(known[_edge(base, v)] + length - known[_edge(v, vertex)]) / 2 for v in others])
            coefficients = (adjugate * h / patch["determinant"]).applyfunc(sp.simplify)
            consistency = sp.simplify((h.T * coefficients)[0] - length)
            if consistency != 0:
                conflicts.append({"kind": "inconsistent_five_point_metric", "anchors": list(anchors),
                                  "vertex": vertex, "squared_norm_discrepancy": consistency})
                continue
            expression = radius_at[vertex] - radius_at[base] - length
            expression += sum(coefficient * (radius_at[base] + known[_edge(base, v)] - radius_at[v])
                              for coefficient, v in zip(coefficients, others))
            index = add_equation(expression, {"kind": "additional_point_affine_radius_relation", "patch": patch_id,
                                             "vertex": vertex, "affine_coefficients": list(coefficients),
                                             "anchors": list(anchors)})
            additional.append({"vertex": vertex, "affine_coefficients": list(coefficients),
                               "exact_metric_consistency": True, "equation_index": index})
        patch["additional_points"] = additional
        metric_patches.append(patch)

    model["patch_equations"] = patch_equations
    model["patch_sources"] = patch_sources
    model["known_distances"] = known
    model["known_distance_sources"] = dict(distance_sources)
    model["reflections"] = reflections
    model["reflection_squares"] = inferred_squares
    model["square_pyramid_patches"] = pyramids
    model["metric_patches"] = metric_patches
    model["patch_metadata"] = {"numerical_positions_used": False, "automorphism_enumeration": symmetry,
                                "cliques_examined": examined, "metric_patches_truncated": truncated,
                                "distance_conflicts": conflicts,
                                "equations_are_necessary_conditions": True,
                                "original_faces_changed": False}
    return model
