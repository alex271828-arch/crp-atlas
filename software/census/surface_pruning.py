"""Exact surface rotation and coupled boundary-link pruning (v1.11)."""

from face_types import TOTAL_CURVATURE, defect

from functools import lru_cache
import boundary_filter as BF


def rotation_faces(faces, info):
    """Return oriented complementary face walks, or None for positive genus.

    At most two partial fans fit a degree-at-most-five final link, since each
    requires one intervening missing face. Their order around a vertex is
    unique up to cyclic shift. Concatenating them therefore gives the unique
    rotation of the existing edges induced by any possible target sphere.
    """
    turn = {}
    for v, (_, parts, _) in info["domains"].items():
        if len(parts) > 2:
            raise ValueError("Rotation uniqueness requires at most two link components")
        nxt = {a: b for a, b, n in info["arcs"][v]}
        order = []
        for start, end, word in parts:
            u = start
            order.append(u)
            for _ in word:
                u = nxt[u]
                if u != start:
                    order.append(u)
        for i, u in enumerate(order):
            turn[(u, v)] = (v, order[(i + 1) % len(order)])
    unseen = set(turn)
    walks = []
    while unseen:
        first = min(unseen)
        e = first
        walk = []
        while e in unseen:
            unseen.remove(e)
            walk.append(e[0])
            e = turn[e]
        if e != first:
            raise ValueError("Invalid rotation permutation")
        walks.append(tuple(walk))
    # Every state is connected by construction.
    if len(info["vertices"]) - len(turn) // 2 + len(walks) != 2:
        return None

    def canon(f):
        return min(f[i:] + f[:i] for i in range(len(f)))

    present = {canon(tuple(f)) for f in faces}
    return tuple(w for w in walks if canon(w) not in present)


@lru_cache(maxsize=50000)
def gap_rows(paths, menu):
    """All joint outside words, in the unique cyclic order of the old fans."""
    rows = set()
    for sig in menu:
        variants = {w[i:] + w[:i] for w in (sig, sig[::-1]) for i in range(len(w))}
        curvature = defect(sig)
        a = len(paths[0])
        for word in variants:
            if word[:a] != paths[0]:
                continue
            if len(paths) == 1:
                if len(word) > a:
                    rows.add(((word[a:],), curvature))
            elif len(paths) == 2:
                b = len(paths[1])
                for k in range(a + 1, len(word) - b):
                    if word[k : k + b] == paths[1]:
                        rows.add(((word[a:k], word[k + b :]), curvature))
            else:
                raise ValueError("At most two fans are supported")
    return tuple(sorted(rows))


@lru_cache(maxsize=50000)
def run_possible(rows):
    """Existence-only forced-run DP; stop after the first supported cycle."""
    n = len(rows)
    if n in BF.ANGLE and all((n, n, 1) in r for r in rows):
        return True
    for cut, options in enumerate(rows):
        seq = rows[cut + 1 :] + rows[:cut]
        for first, out, length in options:
            if length == 1:
                continue
            reachable = {(out, 1)}
            for opts in seq:
                following = set()
                for a, count in reachable:
                    for incoming, b, length in opts:
                        if incoming != a:
                            continue
                        if length == 1:
                            if count + 1 < a:
                                following.add((b, count + 1))
                        else:
                            following.add((b, 1))
                reachable = following
                if not reachable:
                    break
            if any(a == first for a, count in reachable):
                return True
    return False


def refine(faces, info, menu):
    """Joint word consistency across all complementary regions, including pinches.

    Rows couple both gaps at a pinched vertex. Cyclic DP propagates common
    missing-face sizes around each region; row filtering feeds deductions
    back across regions. Curvature is counted once per actual vertex.
    """
    cycles = rotation_faces(faces, info)
    if cycles is None:
        return {"excluded": True, "reason": "rotation_genus"}
    ds = {
        v: gap_rows(tuple(p[2] for p in info["domains"][v][1]), menu)
        for v in info["opened"]
    }
    # A region corner identifies one particular gap, not just a vertex.
    corner_cycles = []
    for cycle in cycles:
        corners = []
        for j, v in enumerate(cycle):
            a, b = cycle[j - 1], cycle[(j + 1) % len(cycle)]
            parts = info["domains"][v][1]
            gap = next(
                (
                    i
                    for i, p in enumerate(parts)
                    if p[1] == a and parts[(i + 1) % len(parts)][0] == b
                ),
                None,
            )
            if gap is None:
                raise ValueError("Region corner is not an unfinished link gap")
            corners.append((v, gap))
        corner_cycles.append(corners)
    while True:
        if any(not rows for rows in ds.values()):
            return {"excluded": True, "reason": "coupled_boundary_words"}
        projected = {
            (v, i): tuple(sorted({(g[i][0], g[i][-1], 0, g[i]) for g, k in rows}))
            for v, rows in ds.items()
            for i in range(len(info["domains"][v][1]))
        }
        supported = BF.supported_domains(corner_cycles, projected)
        allowed = {c: {d[3] for d in rows} for c, rows in supported.items()}
        new = {
            v: tuple(
                r for r in rows if all(g in allowed[v, i] for i, g in enumerate(r[0]))
            )
            for v, rows in ds.items()
        }
        if new == ds:
            break
        ds = new
    debt = sum(defect(d[2][0]) for v, d in info["domains"].items() if d[0])
    debt += sum(min(k for gaps, k in rows) for rows in ds.values())
    if debt > TOTAL_CURVATURE:
        return {"excluded": True, "reason": "coupled_curvature"}
    for cycle in corner_cycles:
        rows = tuple(
            tuple(sorted({(a, b, len(g)) for a, b, k, g in supported[c]}))
            for c in cycle
        )
        if not run_possible(rows):
            return {"excluded": True, "reason": "coupled_face_run"}
    info["gap_rows"] = ds
    info["region_cycles"] = cycles
    info["boundary_debt"] = debt
    return {"excluded": False, "reason": "unknown"}


def corner_allowed(info, a, v, b, n):
    """Use globally supported gap rows before extending a partial new polygon."""
    if v not in info.get("gap_rows", {}):
        return True
    parts = info["domains"][v][1]
    neighbors = {u for x, y, m in info["arcs"][v] for u in (x, y)}
    for gaps, _curvature in info["gap_rows"][v]:
        for i, gap in enumerate(gaps):
            left = parts[i][1]
            right = parts[(i + 1) % len(parts)][0]
            if a in neighbors and b in neighbors:
                if a == left and b == right and gap == (n,):
                    return True
            elif a in neighbors:
                if a == left and len(gap) >= 2 and gap[0] == n:
                    return True
            elif b in neighbors:
                if b == right and len(gap) >= 2 and gap[-1] == n:
                    return True
            elif n in gap[1:-1]:
                return True
    return False


def forced_face(info):
    """A polygon whose complete boundary is forced by supported gap words."""
    rows = info["gap_rows"]
    for cycle in info["region_cycles"]:
        n = len(cycle)
        gapsets = []
        for j, v in enumerate(cycle):
            parts = info["domains"][v][1]
            gap = next(i for i, p in enumerate(parts) if p[1] == cycle[j - 1])
            gapsets.append({r[0][gap] for r in rows[v]})
        single = [all(len(g) == 1 for g in gs) for gs in gapsets]
        if all(single):
            sizes = set.intersection(*({g[0] for g in gs} for gs in gapsets))
            if sizes == {n}:
                return tuple(cycle)
        for i in range(n):
            if single[i]:
                continue
            path = [cycle[i], cycle[(i + 1) % n]]
            j = (i + 1) % n
            sizes = {g[-1] for g in gapsets[i]} & {g[0] for g in gapsets[j]}
            while single[j] and len(path) <= n:
                sizes &= {g[0] for g in gapsets[j]}
                q = (j + 1) % n
                if q == i:
                    break
                path.append(cycle[q])
                sizes &= {g[0] for g in gapsets[q]}
                sizes = {m for m in sizes if m >= len(path)}
                if sizes == {len(path)}:
                    return tuple(path)
                j = q
    return None


def propagate(faces, info, menu):
    """Return the mandatory-face closure, preserving every target completion.

    Every iteration adds a face. The finite curvature bound ensures termination.
    The caller replaces its fragment by this closure; its pivot is unchanged.
    """
    import link_growth as LG

    patch = frozenset(faces)
    current = info
    while True:
        verdict = refine(patch, current, menu)
        if verdict["excluded"]:
            return verdict
        face = forced_face(current)
        if face is None:
            return dict(verdict, faces=patch, info=current)
        if any(BF.intersection_obstruction([old, face]) for old in patch):
            return {"excluded": True, "reason": "forced_face_intersection"}
        patch = patch | {face}
        current = LG.analyze(patch, menu)
        if current is None:
            return {"excluded": True, "reason": "forced_face_link"}
