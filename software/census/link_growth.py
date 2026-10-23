"""Injective, oriented face-complex growth, including disconnected vertex fans.

See EXHAUSTIVENESS.md. All rejection arithmetic here is integer/combinatorial.
No numerical realization, quotient, or simple-boundary assumption is used.
"""

from face_types import TOTAL_CURVATURE, defect

from collections import defaultdict, deque
from functools import lru_cache
import boundary_filter as BF
import surface_pruning as SP

ANGLE = BF.ANGLE  # three-degree units


def rotate(f, v):
    i = f.index(v)
    return f[i:] + f[:i]


def link_components(arcs):
    """arcs=(previous neighbour, next neighbour, face size). None = impossible."""
    nxt, prev = {}, {}
    for a, b, n in arcs:
        if a == b or a in nxt or b in prev:
            return None
        nxt[a] = (b, n)
        prev[b] = a
    unseen = set(nxt) | set(prev)
    paths = []
    for start in sorted(unseen - set(prev)):
        word, u = [], start
        while u in nxt:
            unseen.discard(u)
            u, n = nxt[u]
            word.append(n)
        unseen.discard(u)
        paths.append((start, u, tuple(word)))
    if unseen:
        if paths:
            return None  # closed cycle plus a separate component
        start = min(unseen)
        word, u = [], start
        while u in unseen:
            unseen.remove(u)
            u, n = nxt[u]
            word.append(n)
        if u != start or unseen:
            return None
        return True, ((start, start, tuple(word)),)
    return False, tuple(paths)


@lru_cache(maxsize=50000)
def supported(paths, menu):
    """Pack all directed paths into disjoint runs of one allowed cyclic word.

    Different paths require at least one missing face between them: their
    endpoints are distinct existing edges and may not be identified later.
    """
    total = sum(map(len, paths))
    answer = []
    for sig in menu:
        d = len(sig)
        if total + len(paths) > d:
            continue
        variants = {sig, sig[::-1]}
        good = False
        for word in variants:
            options = []
            for path in paths:
                opts = []
                for k in range(d):
                    if all(word[(k + j) % d] == n for j, n in enumerate(path)):
                        mask = sum(1 << ((k + j) % d) for j in range(len(path)))
                        guard = (
                            mask | (1 << ((k - 1) % d)) | (1 << ((k + len(path)) % d))
                        )
                        opts.append((mask, guard))
                options.append(opts)

            def place(i, used, guarded):
                if i == len(options):
                    return True
                return any(
                    not (mask & guarded)
                    and not (guard & used)
                    and place(i + 1, used | mask, guarded | guard)
                    for mask, guard in options[i]
                )

            if place(0, 0, 0):
                good = True
                break
        if good:
            answer.append(sig)
    return tuple(answer)


def vertex_domain(arcs, menu):
    result = link_components(arcs)
    if result is None:
        return None
    closed, parts = result
    if closed:
        sig = BF.canonical(parts[0][2])
        ds = (sig,) if sig in menu else ()
    else:
        ds = supported(tuple(sorted(p[2] for p in parts)), menu)
    return (closed, parts, ds) if ds else None


def analyze(faces, menu):
    arcs = defaultdict(list)
    directed = {}
    vertices = set()
    for fi, f in enumerate(faces):
        if len(f) not in ANGLE or len(set(f)) != len(f):
            return None
        vertices.update(f)
        for i, v in enumerate(f):
            edge = (v, f[(i + 1) % len(f)])
            if edge in directed:
                return None
            directed[edge] = fi
            arcs[v].append((f[i - 1], f[(i + 1) % len(f)], len(f)))
    domains = {v: vertex_domain(tuple(a), menu) for v, a in arcs.items()}
    if any(d is None for d in domains.values()):
        return None
    debt = sum(min(defect(s) for s in d[2]) for d in domains.values())
    if debt > TOTAL_CURVATURE:
        return None
    opened = {v for v, d in domains.items() if not d[0]}
    pinched = {v for v in opened if len(domains[v][1]) > 1}
    return dict(
        arcs=arcs,
        directed=directed,
        domains=domains,
        vertices=vertices,
        opened=opened,
        pinched=pinched,
        debt=debt,
    )


def choose_pivot(state):
    """Choice affects order, never coverage: every chosen link is fully enumerated."""
    return min(
        state["opened"],
        key=lambda v: (
            min(len(s) - len(state["arcs"][v]) for s in state["domains"][v][2]),
            len(state["domains"][v][2]),
            -len(state["arcs"][v]),
            v,
        ),
    )


def next_faces(faces, pivot, menu, state=None, stats=None):
    """All faces filling one chosen unfinished radial edge at pivot.

    Fresh labels are canonical by order of first appearance. Old closed vertices
    cannot receive a new face. All other old vertices are considered, even when
    they have multiple link components. No partial-face coordinates are used.
    """
    faces = tuple(faces)
    state = state or analyze(faces, menu)
    if state is None or pivot not in state["opened"]:
        return
    # At the end of a directed link path, face v->u is already present;
    # the missing face must use u->v. That face is unique in any target.
    u = min(p[1] for p in state["domains"][pivot][1])
    start = (u, pivot)
    covered = state["directed"]
    base = max(state["vertices"]) + 1
    old = sorted(state["opened"])
    vmax = TOTAL_CURVATURE // min(defect(s) for s in menu)
    corner_cache = {}

    def corner(a, v, b, n):
        key = (a, v, b, n)
        if key not in corner_cache:
            corner_cache[key] = (
                SP.corner_allowed(state, a, v, b, n)
                and vertex_domain(tuple(state["arcs"].get(v, ())) + ((a, b, n),), menu)
                is not None
            )
        return corner_cache[key]

    oldsets = [set(f) for f in faces]

    def intersection_ok(prefix, n):
        # All shared vertices must eventually form at most one common edge.
        for f, vs in zip(faces, oldsets):
            common = [(i, v) for i, v in enumerate(prefix) if v in vs]
            if len(common) > 2:
                return False
            if len(common) == 2:
                (i, a), (j, b) = common
                if (f.index(a) - f.index(b)) % len(f) not in (1, len(f) - 1):
                    return False
                if j - i != 1 and not (i == 0 and j == n - 1):
                    return False
        return True

    def build(prefix, n, nfresh):
        if len(prefix) == n:
            a, b = prefix[-1], prefix[0]
            if (a, b) in covered:
                return
            if not corner(prefix[-2], a, b, n) or not corner(a, b, prefix[1], n):
                return
            if stats is not None:
                stats["faces_generated"] = stats.get("faces_generated", 0) + 1
            yield prefix
            return
        choices = old + (
            [base + nfresh] if len(state["vertices"]) + nfresh < vmax else []
        )
        for w in choices:
            if w in prefix or (prefix[-1], w) in covered:
                continue
            if not corner(prefix[-2], prefix[-1], w, n):
                continue
            extended = prefix + (w,)
            if not intersection_ok(extended, n):
                continue
            yield from build(extended, n, nfresh + (w >= base))

    sizes = sorted({n for sig in state["domains"][pivot][2] for n in sig})
    for n in sizes:
        yield from build(start, n, 0)


def canonical_code(faces, root=None):
    """Exact oriented face-complex code, invariant under relabelling and reflection.

    Enumerate invariantly selected starting face corners, then traverse the dual
    graph. It is connected by construction. The emitted face list reconstructs
    the complex; the root label distinguishes inner states with different pivots.
    No hash is used for equality/deduplication.
    """
    faces = tuple(map(tuple, faces))
    incident = defaultdict(list)
    for f in faces:
        for v in f:
            incident[v].append(len(f))
    signatures = {v: (v == root, tuple(sorted(ns))) for v, ns in incident.items()}
    anchors = []
    for hand in (faces, tuple(tuple(reversed(f)) for f in faces)):
        edge_face = {
            (f[j], f[(j + 1) % len(f)]): i
            for i, f in enumerate(hand)
            for j in range(len(f))
        }
        for i, f in enumerate(hand):
            for k in range(len(f)):
                ff = f[k:] + f[:k]
                anchors.append(
                    ((len(f), tuple(signatures[v] for v in ff)), hand, edge_face, i, ff)
                )
    best_signature = min(a[0] for a in anchors)
    best = None
    for signature, hand, edge_face, first, ff in anchors:
        if signature != best_signature:
            continue
        seen = {first}
        queue = deque([(first, ff)])
        labels = {}
        encoded = []
        while queue:
            i, f = queue.popleft()
            for v in f:
                if v not in labels:
                    labels[v] = len(labels)
            encoded.append(tuple(labels[v] for v in f))
            for j, a in enumerate(f):
                b = f[(j + 1) % len(f)]
                other = edge_face.get((b, a))
                if other is not None and other not in seen:
                    seen.add(other)
                    queue.append((other, rotate(hand[other], b)))
        if len(seen) != len(faces):
            raise ValueError("Canonicalization requires an edge-connected face complex")
        code = (labels[root] if root is not None else -1, tuple(encoded))
        if best is None or code < best:
            best = code
    return best
