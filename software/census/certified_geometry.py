"""Conservative convex geometry in fixed-point, outward-rounded integer intervals.

No floating-point geometric arithmetic, sampled tables, acos, epsilon cutoffs,
or midpoint choices of geometric branches. Time budgets only return UNKNOWN.
Angles are represented by cosines on [-1,1]; convexity fixes nonnegative sine.
"""

from collections import defaultdict, deque
from functools import lru_cache
from math import isqrt
import time
from geometry_plan import plan  # purely combinatorial choice of free edge variables

S = 1 << 52


class Impossible(Exception):
    pass


class I:
    """The closed interval [a/S, b/S]; endpoints and scalar operands are integers."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=None):
        self.a = a
        self.b = a if b is None else b

    @staticmethod
    def coerce(x):
        return x if isinstance(x, I) else I(x * S)

    def __add__(self, x):
        x = I.coerce(x)
        return I(self.a + x.a, self.b + x.b)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.b, -self.a)

    def __sub__(self, x):
        return self + -I.coerce(x)

    def __rsub__(self, x):
        return I.coerce(x) + -self

    def __mul__(self, x):
        x = I.coerce(x)
        p = (self.a * x.a, self.a * x.b, self.b * x.a, self.b * x.b)
        return I(min(p) // S, -(-max(p) // S))

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = I.coerce(x)
        if x.a <= 0 <= x.b:
            raise ZeroDivisionError
        q = [(a * S, b) for a in (self.a, self.b) for b in (x.a, x.b)]
        return I(min(a // b for a, b in q), max(-((-a) // b) for a, b in q))

    def width(self):
        return self.b - self.a


FULL = I(-S, S)


def intersect(a, b):
    lo = max(a.a, b.a)
    hi = min(a.b, b.b)
    if lo > hi:
        raise Impossible
    return I(lo, hi)


def cosine(x):
    return intersect(x, FULL)


def square(x):
    a = 0 if x.a <= 0 <= x.b else min(x.a * x.a, x.b * x.b)
    b = max(x.a * x.a, x.b * x.b)
    return I(a // S, -(-b // S))


def sqrt(x):
    if x.b < 0:
        raise Impossible
    lo = isqrt(max(0, x.a) * S)
    hi = isqrt(x.b * S)
    return I(lo, hi + (hi * hi != x.b * S))


@lru_cache(maxsize=65536)
def _sine_bounds(lo, hi):
    """Cache by exact endpoints; no quantization or approximation is introduced."""
    return sqrt(1 - square(cosine(I(lo, hi))))


def sine(c):
    return _sine_bounds(c.a, c.b)


SQRT5 = sqrt(I(5 * S))
SQRT2 = sqrt(I(2 * S))
SIDE = {
    3: I(S // 2),
    4: I(0),
    5: (1 - SQRT5) / 4,
    6: I(-S // 2),
    8: -SQRT2 / 2,
    10: -(1 + SQRT5) / 4,
}
TURN = {
    3: I(-S // 2),
    4: I(0),
    5: (SQRT5 - 1) / 4,
    6: I(S // 2),
    8: SQRT2 / 2,
    10: (1 + SQRT5) / 4,
}
SIN_SIDE = {n: sine(c) for n, c in SIDE.items()}
SIN_TURN = {n: sine(c) for n, c in TURN.items()}

EQUAL_SIDE = {id(c): (c / SIN_SIDE[n], (2 * square(c) - 1).a) for n, c in SIDE.items()}


def tri(opp, a, b):
    """Cosine of the angle opposite side opp; arguments are side cosines.

    Spherical cosine law: cos(A) = (cos(c)-cos(a)cos(b))/(sin(a)sin(b)).
    A singular denominator returns an outer bound, never a guessed branch.
    """
    # Cancel a removable singularity only for IDENTICAL algebraic side
    # constants, never for merely overlapping or equal-looking intervals.
    if opp is a and id(a) in EQUAL_SIDE:
        side = a
        diag = b
    elif opp is b and id(b) in EQUAL_SIDE:
        side = b
        diag = a
    else:
        side = None
    if side is not None:
        cot, lo = EQUAL_SIDE[id(side)]
        if cot.a == cot.b == 0:
            return I(0)
        diag = intersect(diag, I(lo, S))
        return cosine(cot * sqrt((1 - diag) / (1 + diag)))
    den = sine(a) * sine(b)
    numerator = opp - a * b
    # numerator = denominator * cos(angle), including singular triangles.
    if numerator.a > den.b or numerator.b < -den.b:
        raise Impossible
    # Division may be singular, but its undivided equation was still checked.
    if den.a <= 0:
        return FULL
    return cosine(numerator / den)


def angle_sum(a, b):
    """Enclose cos(A+B) monotonically on the convex domain A+B <= pi."""
    if a.b + b.b < 0:
        raise Impossible

    def endpoint(x, y):
        x, y = I(x), I(y)
        return x * y - sine(x) * sine(y)

    lo = endpoint(a.a, b.a).a if a.a + b.a >= 0 else -S
    hi = endpoint(a.b, b.b).b
    return cosine(I(lo, hi))


def angle_difference(a, b):
    """Enclose cos(A-B) monotonically on the domain A >= B."""
    if a.a > b.b:
        raise Impossible

    def endpoint(x, y):
        x, y = I(x), I(y)
        return x * y + sine(x) * sine(y)

    lo = endpoint(a.a, b.b).a
    hi = endpoint(a.b, b.a).b if a.b <= b.a else S
    return cosine(I(lo, hi))


@lru_cache(None)
def pins(word):
    s = [SIDE[n] for n in word]
    return tuple(tri(s[(k + 2) % 3], s[k], s[(k + 1) % 3]) for k in range(3))


def _quad_raw(word, k, c0):
    """Split a convex spherical quadrilateral along the diagonal opposite angle k.

    c0 is cos(angle k). Triangle cosine laws give the opposite angle;
    the two pairs of subangles add to the remaining boundary angles.
    """
    a, b, c, d = [SIDE[word[(k + j) % 4]] for j in range(4)]
    diag = cosine(a * b + sine(a) * sine(b) * c0)
    tb = tri(a, b, diag)
    ta = tri(b, a, diag)
    tc = tri(d, c, diag)
    td = tri(c, d, diag)
    result = (c0, angle_sum(tb, tc), tri(diag, c, d), angle_sum(ta, td))
    return {((k + j) % 4): r for j, r in enumerate(result)}


@lru_cache(maxsize=65536)
def _quad_cached(word, lo, hi):
    try:
        r = _quad_raw(word, 0, I(lo, hi))
        return tuple(r[k] for k in range(4))
    except Impossible:
        return None


def quad(word, k, c0):
    word = tuple(word)
    r = _quad_cached(word[k:] + word[:k], c0.a, c0.b)
    if r is None:
        raise Impossible
    return {(k + j) % 4: c for j, c in enumerate(r)}


@lru_cache(None)
def quad_domain(rotated_word):
    """Outer hull from a complete interval partition, never point sampling."""
    good = []
    for i in range(128):
        a = -S + 2 * S * i // 128
        b = -S + 2 * S * (i + 1) // 128
        try:
            quad(rotated_word, 0, I(a, b))
            good.append((a, b))
        except Impossible:
            pass
    if not good:
        raise Impossible
    return I(min(a for a, b in good), max(b for a, b in good))


def _pent_raw(sizes, a0, a1):
    """Triangulate a convex spherical pentagon from adjacent angle cosines.

    d joins rays 4 and 1; e joins rays 4 and 2. Convexity fixes the
    nonnegative residual at ray 1 and the sums at rays 2 and 4.
    """
    s0, s1, s2, s3, s4 = [SIDE[n] for n in sizes]
    d = cosine(s0 * s4 + sine(s0) * sine(s4) * a0)
    beta = tri(s4, s0, d)
    g4a = tri(s0, s4, d)
    # Residual angle a1-beta must be nonnegative.
    residual = angle_difference(a1, beta)
    e = cosine(s1 * d + sine(s1) * sine(d) * residual)
    q2a = tri(d, s1, e)
    q4b = tri(s1, d, e)
    q2b = tri(s3, s2, e)
    q4c = tri(s2, s3, e)
    q3 = tri(e, s2, s3)
    return (a0, a1, angle_sum(q2a, q2b), q3, angle_sum(angle_sum(g4a, q4b), q4c))


@lru_cache(maxsize=65536)
def _pent_cached(word, a, b, c, d):
    try:
        return _pent_raw(word, I(a, b), I(c, d))
    except Impossible:
        return None


def pent(word, x, y):
    r = _pent_cached(tuple(word), x.a, x.b, y.a, y.b)
    if r is None:
        raise Impossible
    return r


def vertex_words(info):
    out = {}
    for v, (closed, parts, ds) in info["domains"].items():
        if not closed:
            continue
        nxt = {a: (b, n) for a, b, n in info["arcs"][v]}
        start = parts[0][0]
        u = start
        word = []
        edges = []
        while True:
            u, n = nxt[u]
            word.append(n)
            edges.append(frozenset((v, u)))
            if u == start:
                break
        out[v] = (tuple(word), edges)
    return out


def refine_figures(vw, angles):
    angles = dict(angles)

    def cells(x, m):
        return [
            I(x.a + x.width() * i // m, x.a + x.width() * (i + 1) // m)
            for i in range(m)
        ]

    for w, es in vw.values():
        results = []
        if len(w) == 4:
            k = min(range(4), key=lambda k: angles[es[k]].width())
            x = angles[es[k]]
            for cell in cells(x, 12 if x.width() > S // 100 else 1):
                try:
                    results.append(
                        {
                            j: intersect(c, angles[es[j]])
                            for j, c in quad(w, k, cell).items()
                        }
                    )
                except Impossible:
                    pass
        elif len(w) == 5:
            k0 = min(
                range(5),
                key=lambda k: angles[es[k]].width() + angles[es[(k + 1) % 5]].width(),
            )
            k1 = (k0 + 1) % 5
            rot = tuple(w[(k1 + j) % 5] for j in range(5))
            x, y = angles[es[k0]], angles[es[k1]]
            for a in cells(x, 5 if x.width() > S // 100 else 1):
                for b in cells(y, 5 if y.width() > S // 100 else 1):
                    try:
                        r = pent(rot, a, b)
                        results.append(
                            {
                                (k1 + j - 1)
                                % 5: intersect(c, angles[es[(k1 + j - 1) % 5]])
                                for j, c in enumerate(r)
                            }
                        )
                    except Impossible:
                        pass
        else:
            continue
        if not results:
            raise Impossible
        for k, e in enumerate(es):
            angles[e] = I(min(r[k].a for r in results), max(r[k].b for r in results))
    return angles


def evaluate(vw, order, box):
    angles = {e: FULL for word, es in vw.values() for e in es}
    for w, es in vw.values():
        if len(w) == 4:
            for k, e in enumerate(es):
                angles[e] = intersect(angles[e], quad_domain(w[k:] + w[:k]))
    for e, c in box.items():
        angles[e] = intersect(angles.get(e, FULL), c)
    for v, ks in order:
        word, es = vw[v]
        if len(word) == 3:
            result = dict(enumerate(pins(word)))
        elif len(word) == 4:
            k = min(range(4), key=lambda k: angles[es[k]].width())
            result = quad(word, k, angles[es[k]])
        else:
            k0, k1 = ks
            rot = tuple(word[(k1 + j) % 5] for j in range(5))
            r = pent(rot, angles[es[k0]], angles[es[k1]])
            result = {(k1 + j - 1) % 5: c for j, c in enumerate(r)}
        for k, c in result.items():
            angles[es[k]] = intersect(angles[es[k]], c)
    return refine_figures(vw, angles)


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(c, a):
    return [c * x for x in a]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), I(0))


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def unit_components(a):
    return [cosine(x) for x in a]


def convex_check(faces, angles):
    """Develop oriented regular polygons; convex folds have a fixed sign.

    The initial face is in z=0 with outward normal +z. A reflection of the
    target supplies this orientation if needed. No midpoint sign is chosen.
    """
    faces = tuple(faces)
    ef = defaultdict(list)
    for fi, f in enumerate(faces):
        for j, v in enumerate(f):
            ef[frozenset((v, f[j - 1]))].append(fi)
    X = {}
    N = {0: [I(0), I(0), I(S)]}

    def put(v, p):
        X[v] = [intersect(a, b) for a, b in zip(X[v], p)] if v in X else p
        return X[v]

    def polygon(f, start, direction, normal):
        pos = start
        put(f[0], pos)
        n = len(f)
        for v in f[1:]:
            pos = put(v, add(pos, direction))
            direction = unit_components(
                add(
                    scale(TURN[n], direction),
                    scale(SIN_TURN[n], cross(normal, direction)),
                )
            )

    polygon(faces[0], [I(0), I(0), I(0)], [I(S), I(0), I(0)], N[0])
    placed = {0}
    queue = deque([0])
    while queue:
        fi = queue.popleft()
        f = faces[fi]
        nf = N[fi]
        for j, a in enumerate(f):
            b = f[(j + 1) % len(f)]
            e = frozenset((a, b))
            for gj in ef[e]:
                if gj in placed:
                    continue
                c = angles.get(e, FULL)
                s = sine(c)
                u = unit_components(sub(X[b], X[a]))
                inside = unit_components(cross(nf, u))
                ng = unit_components(sub(scale(-s, inside), scale(c, nf)))
                g = faces[gj]
                k = g.index(b)
                g = g[k:] + g[:k]
                if g[1] != a:
                    raise ValueError("Inconsistent face orientation")
                polygon(g, X[b], scale(-1, u), ng)
                N[gj] = ng
                placed.add(gj)
                queue.append(gj)
    for fi, f in enumerate(faces):
        for v, p in X.items():
            if v not in f and dot(N[fi], sub(p, X[f[0]])).a > 0:
                raise Impossible


def solve(faces, info, parent=None, time_cap=0.1, max_boxes=200000, cover_limit=128):
    """Empty cover certifies exclusion; every timeout retains the full remaining cover."""
    import rim_geometry as RG

    start = time.monotonic()
    faces = tuple(sorted(map(tuple, faces)))
    rim = RG.fans(info)
    vw = vertex_words(info)
    if not vw:
        return False, None, 0
    order, free = plan(vw)
    inherited = parent is not None and set(parent["faces"]).issubset(faces)
    serial = parent["boxes"] if inherited else [{}]
    boxes = []
    for old in serial:
        box = {e: FULL for e in free}
        box.update({e: I(a, b) for e, (a, b) in old.items()})
        boxes.append(box)
    count = 0
    unresolved = []

    def finish():
        all_boxes = unresolved + boxes
        saved = [{e: (x.a, x.b) for e, x in b.items()} for b in all_boxes]
        if len(saved) > cover_limit:
            keys = set.intersection(*(set(b) for b in saved))
            saved = [
                {
                    e: (min(b[e][0] for b in saved), max(b[e][1] for b in saved))
                    for e in keys
                }
            ]
        return not saved, {"faces": faces, "boxes": saved}, count

    while boxes:
        if time.monotonic() - start >= time_cap or count >= max_boxes:
            return finish()
        bx = boxes.pop()
        count += 1
        try:
            for _ in range(4):
                result = evaluate(vw, order, bx)
                RG.check(rim, result)
                shrink = sum(bx[e].width() - result[e].width() for e in bx)
                bx = result
                if shrink < S // 1000000:
                    break
            # Check the entire coordinate enclosure, without selecting a point.
            convex_check(faces, result)
        except Impossible:
            continue
        if free:
            e = max(free, key=lambda e: bx[e].width())
            if bx[e].width() > S // 1000000:
                mid = (bx[e].a + bx[e].b) // 2
                left = dict(bx)
                right = dict(bx)
                left[e] = I(bx[e].a, mid)
                right[e] = I(mid, bx[e].b)
                boxes.extend((left, right))
                continue
        unresolved.append(result)
    return finish()
