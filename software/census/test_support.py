"""Regression-only reference routines; not imported by the census."""

from face_types import ANGLE
from boundary_filter import canonical


def vertices_of(faces):
    return {v for f in faces for v in f}


def link_signature(faces, vertex):
    arcs = {
        f[f.index(vertex) - 1]: (f[(f.index(vertex) + 1) % len(f)], len(f))
        for f in faces
        if vertex in f
    }
    start = min(arcs)
    word = []
    v = start
    while True:
        v, n = arcs[v]
        word.append(n)
        if v == start:
            return canonical(tuple(word))


def cycle_minimum(cycle, ds, run_limits=False):
    if not run_limits:
        best = 10**9
        for first in ANGLE:
            cost = {first: 0}
            for v in cycle:
                nxt = {}
                for a, b, k, _ in ds[v]:
                    if a in cost:
                        nxt[b] = min(nxt.get(b, 10**9), cost[a] + k)
                cost = nxt
            best = min(best, cost.get(first, 10**9))
        return best
    # Across a vertex with one missing face, the two boundary edges belong to
    # the SAME future polygon. Track its accumulated boundary-edge count.
    # Any proper run in an m-gon has at most m-1 edges. A run covering the
    # entire boundary is possible only when that boundary itself is the m-gon.
    best = 10**9
    n = len(cycle)
    for m in ANGLE:
        if n == m:
            choices = [[k for a, b, k, rest in ds[v] if rest == (m,)] for v in cycle]
            if all(choices):
                best = min(best, sum(min(c) for c in choices))
    # Cut at a corner where at least two faces are missing. Such a corner
    # separates two distinct new faces and makes run counts unambiguous.
    for cut in range(n):
        seq = cycle[cut:] + cycle[:cut]
        v = seq[0]
        for first, out, k0, rest in ds[v]:
            if len(rest) == 1:
                continue
            cost = {(out, 1): k0}
            for w in seq[1:]:
                nxt = {}
                for (a, count), cost0 in cost.items():
                    for incoming, b, k, r in ds[w]:
                        if incoming != a:
                            continue
                        if len(r) == 1:
                            if count + 1 >= a:
                                continue
                            key = (b, count + 1)
                        else:
                            key = (b, 1)
                        nxt[key] = min(nxt.get(key, 10**9), cost0 + k)
                cost = nxt
            best = min(
                best,
                min((c for (a, count), c in cost.items() if a == first), default=10**9),
            )
    return best
