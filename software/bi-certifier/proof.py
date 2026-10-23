"""Shared exact proof checks and the square contraction certificate.

Only the producer uses approximate arithmetic. A certificate supplies rational
numbers; every proof inequality below is evaluated with outward intervals.
"""

from fractions import Fraction
import exact_intervals as iv


class CertificateError(AssertionError):
    """A sufficient proof condition was not established; never an exclusion."""


def require(condition, message="Certificate condition failed"):
    if not condition:
        raise CertificateError(message)


def rational(value):
    require(
        type(value) in (str, int, Fraction), "Expected an exact rational, not a float or boolean"
    )
    return Fraction(value)


def precision(bits):
    require(type(bits) is int and 80 <= bits <= 4096, "bits outside 80..4096")
    iv.BITS = bits
    iv.SCALE = 1 << bits


def root_box(cert, dimension):
    require(
        isinstance(cert["center"], list) and len(cert["center"]) == dimension and dimension > 0,
        "Invalid root dimension",
    )
    center = [iv.I(rational(x)) for x in cert["center"]]
    radius = iv.I(rational(cert["radius"]))
    require(radius.lo > 0, "Root radius must be positive")
    box = [iv.I(c.lo - radius.hi, c.hi + radius.hi, raw=True) for c in center]
    return center, box


def contraction_image(cert, center, box, residual, jacobian):
    """Enclose c-YF(c)+(I-YJ(X))(X-c), and prove strict contraction.

    Jacobian rows are sparse maps from column to outward interval. Exact input
    constants need only be enclosed: the bounds cover the fixed true constants.
    The same argument applies to a continuous finite piecewise-C1 map, provided
    the caller bounds every piece and proves agreement at transitions.
    """
    n = len(center)
    require(len(residual) == len(jacobian) == n, "Equation dimension mismatch")
    data = cert["preconditioner"]
    require(
        isinstance(data, list)
        and len(data) == n
        and all(isinstance(row, list) and len(row) == n for row in data),
        "Invalid preconditioner",
    )
    inverse = [[iv.I(rational(x)) for x in row] for row in data]
    delta = [x - c for x, c in zip(box, center)]
    # Transpose the sparse Jacobian once; absent entries are exact zeros.
    columns = [[] for _ in range(n)]
    for k, derivative in enumerate(jacobian):
        for j, value in derivative.items():
            columns[j].append((k, value))
    image = []
    contraction = 0
    for i, inverse_row in enumerate(inverse):
        row = [
            iv.I(int(i == j)) - iv.sum_products((inverse_row[k], value) for k, value in column)
            for j, column in enumerate(columns)
        ]
        contraction = max(contraction, sum(x.abs_upper() for x in row))
        require(contraction < iv.SCALE, "Contraction bound is not below one")
        value = center[i] - iv.sum_products(zip(inverse_row, residual))
        value += iv.sum_products(zip(row, delta))
        require(value.strictly_inside(box[i]), "Krawczyk inclusion failed")
        image.append(value)
    return image, str(Fraction(contraction, iv.SCALE))
