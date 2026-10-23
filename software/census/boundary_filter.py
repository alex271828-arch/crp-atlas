"""Exact necessary conditions for injective regular-face patch completion.

No geometric tolerances; angle units are three degrees. UNKNOWN is not existence.
A boundary-cycle letter is the size of the missing face across an edge.
Only cyclic support propagation is used; old single-boundary filters are gone.
"""

from itertools import product, combinations

from face_types import ANGLE, FULL_TURN


def canonical(w):
    return min(
        t[i:] + t[:i] for t in (tuple(w), tuple(reversed(w))) for i in range(len(w))
    )


WORDS = sorted(
    {
        canonical(w)
        for d in (3, 4, 5)
        for w in product(ANGLE, repeat=d)
        if sum(ANGLE[x] for x in w) < FULL_TURN
    }
)
# Pass indices belong to the caller's ordering, not to this lexical alphabet.


def intersection_obstruction(faces):
    for i, j in combinations(range(len(faces)), 2):
        a, b = faces[i], faces[j]
        common = set(a) & set(b)
        if len(common) > 2:
            return {"reason": "facets_share_three_vertices", "faces": [i, j]}
        if len(common) == 2:
            x, y = common
            if any(
                (f.index(x) - f.index(y)) % len(f) not in (1, len(f) - 1)
                for f in (a, b)
            ):
                return {"reason": "facets_share_nonedge_pair", "faces": [i, j]}
    return None


def supported_domains(cycles, ds):
    """Remove rows lacking support in the cyclic missing-face-size equations."""
    ds = dict(ds)
    for cyc in cycles:
        support = {v: set() for v in cyc}
        n = len(cyc)
        for first in ANGLE:
            forward = [{first}]
            for v in cyc:
                forward.append({b for a, b, k, r in ds[v] if a in forward[-1]})
            backward = [set() for _ in range(n + 1)]
            backward[n] = {first}
            for i in range(n - 1, -1, -1):
                backward[i] = {a for a, b, k, r in ds[cyc[i]] if b in backward[i + 1]}
            for i, v in enumerate(cyc):
                support[v].update(
                    d for d in ds[v] if d[0] in forward[i] and d[1] in backward[i + 1]
                )
        for v in cyc:
            ds[v] = tuple(sorted(support[v]))
    return ds
