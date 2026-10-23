"""Exact outward dyadic interval arithmetic. No machine floats in proof operations.
Endpoints are integers divided by 2**BITS. Division uses integer floor/ceiling;
sqrt uses math.isqrt. The checker accepts only rational inputs.
"""

from fractions import Fraction
from math import isqrt

BITS = 160
SCALE = 1 << BITS


def ceildiv(a, b):
    return -((-a) // b)


class DomainError(ArithmeticError):
    pass


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo=0, hi=None, *, raw=False):
        if isinstance(lo, I) and hi is None:
            self.lo, self.hi = lo.lo, lo.hi
            return
        if raw:
            if type(lo) is not int or (hi is not None and type(hi) is not int):
                raise TypeError("Raw endpoints must be integers")
            self.lo, self.hi = lo, lo if hi is None else hi
        elif type(lo) is int and (hi is None or type(hi) is int):
            # Integer endpoints are already exact: no Fraction construction or rounding.
            self.lo = lo * SCALE
            self.hi = (lo if hi is None else hi) * SCALE
        else:
            if isinstance(lo, float) or isinstance(hi, float):
                raise TypeError("Proof inputs must be exact rational values, not float")
            a, b = Fraction(lo), Fraction(lo if hi is None else hi)
            if a > b:
                raise ValueError("reversed interval")
            self.lo = (a.numerator * SCALE) // a.denominator
            self.hi = ceildiv(b.numerator * SCALE, b.denominator)
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    def __add__(self, x):
        if not isinstance(x, I):
            x = I(x)
        return I(self.lo + x.lo, self.hi + x.hi, raw=True)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo, raw=True)

    def __sub__(self, x):
        return self + (-I(x))

    def __rsub__(self, x):
        return I(x) + (-self)

    def __mul__(self, x):
        if not isinstance(x, I):
            x = I(x)
        z = (self.lo * x.lo, self.lo * x.hi, self.hi * x.lo, self.hi * x.hi)
        return I(min(z) // SCALE, ceildiv(max(z), SCALE), raw=True)

    __rmul__ = __mul__

    def __truediv__(self, x):
        if not isinstance(x, I):
            x = I(x)
        if x.lo <= 0 <= x.hi:
            raise DomainError("division interval contains zero")
        if x.hi < 0:
            return (-self) / (-x)
        low = min((a * SCALE) // b for a in (self.lo, self.hi) for b in (x.lo, x.hi))
        high = max(ceildiv(a * SCALE, b) for a in (self.lo, self.hi) for b in (x.lo, x.hi))
        return I(low, high, raw=True)

    def __rtruediv__(self, x):
        return I(x) / self

    def square(self):
        a = 0 if self.lo <= 0 <= self.hi else min(self.lo * self.lo, self.hi * self.hi)
        b = max(self.lo * self.lo, self.hi * self.hi)
        return I(a // SCALE, ceildiv(b, SCALE), raw=True)

    def sqrt(self):
        if self.lo < 0:
            raise DomainError("sqrt interval is not nonnegative")
        a = isqrt(self.lo * SCALE)
        n = self.hi * SCALE
        b = isqrt(n)
        if b * b < n:
            b += 1
        return I(a, b, raw=True)

    def abs_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def midpoint(self):
        return I((self.lo + self.hi) // 2, raw=True)

    def strictly_inside(self, x):
        return self.lo > x.lo and self.hi < x.hi

    def contains_zero(self):
        return self.lo <= 0 <= self.hi

    def intersect(self, x):
        lo, hi = max(self.lo, x.lo), min(self.hi, x.hi)
        if lo > hi:
            raise DomainError("empty intersection")
        return I(lo, hi, raw=True)

    def exact_endpoints(self):
        return [str(Fraction(self.lo, SCALE)), str(Fraction(self.hi, SCALE))]

    def display(self):
        return [float(Fraction(self.lo, SCALE)), float(Fraction(self.hi, SCALE))]

    def __repr__(self):
        return f"I{self.display()}"


def sum_products(pairs):
    """Enclose a dot product with one outward rounding, using integer sums.

    Each endpoint product has denominator SCALE**2. Summing their minima
    and maxima before division is exact and avoids rounding every term.
    """
    low = high = 0
    for a, b in pairs:
        products = (a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi)
        low += min(products)
        high += max(products)
    return I(low // SCALE, ceildiv(high, SCALE), raw=True)
