"""Conservative endpoint-distance tests for incomplete spherical vertex links."""

from collections import Counter
from functools import lru_cache
import certified_geometry as C
import boundary_filter as B

from face_types import HALF_TURN

# Exact algebraic cos/sin(3 degrees) = cos/sin(18 - 15 degrees).
c36 = (1 + C.SQRT5) / 4
c18 = C.sqrt((1 + c36) / 2)
s18 = C.sqrt((1 - c36) / 2)
c30 = C.sqrt(C.I(3 * C.S)) / 2
c15 = C.sqrt((1 + c30) / 2)
s15 = C.sqrt((1 - c30) / 2)
c3 = c18 * c15 + s18 * s15
s3 = s18 * c15 - c18 * s15
COS = [C.I(C.S)]
s = C.I(0)
c = COS[0]
for _ in range(HALF_TURN):
    c, s = C.cosine(c * c3 - s * s3), C.cosine(s * c3 + c * s3)
    COS.append(c)
COS[0] = C.I(C.S)
COS[HALF_TURN // 2] = C.I(0)
COS[HALF_TURN] = C.I(-C.S)


@lru_cache(None)
def ranges(path, words):
    out = []
    for w in words:
        remaining = Counter(w) - Counter(path)
        angles = [B.ANGLE[n] for n, count in remaining.items() for _ in range(count)]
        if not angles:
            continue
        total = sum(angles)
        lo = max(0, 2 * max(angles) - total)
        hi = min(HALF_TURN, total)
        if lo <= hi:
            out.append(C.I(COS[hi].a, COS[lo].b))
    return tuple(out)


def fans(info):
    out = []
    for v in info["opened"]:
        closed, parts, words = info["domains"][v]
        nxt = {a: b for a, b, n in info["arcs"][v]}
        for start, end, path in parts:
            if len(path) < 2:
                continue
            u = start
            edges = []
            for _ in path[:-1]:
                u = nxt[u]
                edges.append(frozenset((v, u)))
            out.append((path, edges, ranges(path, words)))
    return out


def chord(path, angles):
    c = C.SIDE[path[0]]
    previous = None
    for k in range(1, len(path)):
        a = angles[k - 1]
        if k > 1:
            beta = C.tri(previous, C.SIDE[path[k - 1]], c)
            a = C.angle_difference(a, beta)
        previous, c = c, C.cosine(
            c * C.SIDE[path[k]] + C.sine(c) * C.SIN_SIDE[path[k]] * a
        )
    return c


def check(fans, angles):
    for path, edges, allowed in fans:
        if any(e not in angles for e in edges):
            continue
        c = chord(path, [angles[e] for e in edges])
        if not any(c.b >= x.a and c.a <= x.b for x in allowed):
            raise C.Impossible
