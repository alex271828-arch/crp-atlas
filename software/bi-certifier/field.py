"""Exact Q(sqrt(2),sqrt(3),sqrt(5)); proposed identities are never inferred by tolerance."""

from fractions import Fraction as Q
from functools import lru_cache
from exact_intervals import I

RAD = (1, 2, 3, 6, 5, 10, 15, 30)


class K:
    def __init__(self, x=0):
        self.a = (
            x.a
            if isinstance(x, K)
            else tuple(map(Q, x)) if isinstance(x, (list, tuple)) else (Q(x),) + (Q(0),) * 7
        )
        if len(self.a) != 8:
            raise ValueError("Field coefficients must have length eight")

    def __add__(self, b):
        return K([x + y for x, y in zip(self.a, K(b).a)])

    __radd__ = __add__

    def __neg__(self):
        return K([-x for x in self.a])

    def __sub__(self, b):
        return self + -K(b)

    def __rsub__(self, b):
        return K(b) + -self

    def __mul__(self, b):
        out = [Q(0)] * 8
        for i, x in enumerate(self.a):
            if x:
                for j, y in enumerate(K(b).a):
                    if y:
                        out[i ^ j] += x * y * RAD[i & j]
        return K(out)

    __rmul__ = __mul__

    def inv(self):
        if not any(self.a[1:]):
            if not self.a[0]:
                raise ZeroDivisionError("Noninvertible field element")
            return K(1 / self.a[0])
        # Eliminate each quadratic generator by conjugation.
        a = self
        numerator = K(1)
        for bit in (1, 2, 4):
            conjugate = K([-v if i & bit else v for i, v in enumerate(a.a)])
            numerator = numerator * conjugate
            a = a * conjugate
        if any(a.a[1:]) or not a.a[0]:
            raise ZeroDivisionError("Noninvertible field element")
        return K([v / a.a[0] for v in numerator.a])

    def __truediv__(self, b):
        return self * K(b).inv()

    def __eq__(self, b):
        return self.a == K(b).a

    def interval(self):
        return sum((I(v) * I(RAD[i]).sqrt() for i, v in enumerate(self.a) if v), I(0))

    def pack(self):
        return list(map(str, self.a))


def radical(d):
    a = [0] * 8
    a[RAD.index(d)] = 1
    return K(a)


@lru_cache(maxsize=64)
def chord_squared(m, d):
    if m not in (3, 4, 5, 6, 8, 10):
        raise ValueError("Unsupported polygon size")
    d = min(d, m - d)
    if d == 0:
        return K(0)
    t = {
        3: K(-1),
        4: K(0),
        5: (radical(5) - 1) / 2,
        6: K(1),
        8: radical(2),
        10: (radical(5) + 1) / 2,
    }[m]
    c0, c1 = K(2), t
    for _ in range(1, d):
        c0, c1 = c1, t * c1 - c0
    return (2 - c1) / (2 - t)
