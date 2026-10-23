"""Exact certificates for an interior-center radial realization.

No floating-point comparison is a proof step here.  The base field is given by
an isolated real primitive element.  Rational interval arithmetic proves signs
and winding; small multiquadratic algebras prove the fan identities exactly.

The geometric certificate is: a spherical cell decomposition with positive
face cones and one turn at every vertex, regular flat faces, and convex hinges.
Developing the cones then gives a local isometry of the sphere, hence a global
one (the input is checked to be a connected combinatorial sphere).  The radial
surface is consequently embedded.  Convex hinges make it a convex polyhedron.
Coplanar adjacent faces are deliberately allowed.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import isqrt

import sympy as sp


class ProofFailure(Exception):
    """A proof obligation could not be established."""


def _q(value):
    value = sp.Rational(value)
    return Fraction(int(value.p), int(value.q))


def _sq(value):
    return sp.Rational(value.numerator, value.denominator)


class Interval:
    """Closed rational intervals with outward-rounded square roots."""

    def __init__(self, lo, hi=None):
        self.lo = Fraction(lo)
        self.hi = Fraction(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ProofFailure("Invalid interval")

    @staticmethod
    def cast(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __add__(self, other):
        other = self.cast(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        other = self.cast(other)
        vals = [self.lo * other.lo, self.lo * other.hi,
                self.hi * other.lo, self.hi * other.hi]
        return Interval(min(vals), max(vals))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ProofFailure("An interval denominator contains zero")
        return Interval(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * self.cast(other).reciprocal()

    def __pow__(self, n):
        if n < 0:
            return self.reciprocal() ** (-n)
        if n == 0:
            return Interval(1)
        if n % 2:
            return Interval(self.lo ** n, self.hi ** n)
        low = 0 if self.lo <= 0 <= self.hi else min(self.lo ** n, self.hi ** n)
        return Interval(low, max(self.lo ** n, self.hi ** n))

    def sqrt(self, digits=36):
        if self.lo < 0:
            raise ProofFailure("A square-root interval is not nonnegative")
        scale = 10 ** digits

        def lower(x):
            return Fraction(isqrt(x.numerator * scale * scale // x.denominator), scale)

        low, high = lower(self.lo), lower(self.hi)
        if high * high != self.hi:
            high += Fraction(1, scale)
        return Interval(low, high)

    def strings(self, digits=18):
        # Keep proof files compact while preserving outward rounding.
        scale = 10 ** digits
        low = self.lo.numerator * scale // self.lo.denominator
        high = -((-self.hi.numerator * scale) // self.hi.denominator)
        return [str(Fraction(low, scale)), str(Fraction(high, scale))]


def _atan_small(x, terms=32):
    """Alternating-series rational enclosure, for 0 <= x <= 1/2."""
    if not 0 <= x <= Fraction(1, 2):
        raise ProofFailure("Internal arctangent argument outside its range")
    total = Fraction(0)
    power = x
    for k in range(terms):
        total += (-1 if k % 2 else 1) * power / (2 * k + 1)
        power *= x * x
    following = (-1 if terms % 2 else 1) * power / (2 * terms + 1)
    return Interval(min(total, total + following), max(total, total + following))


PI = 16 * _atan_small(Fraction(1, 5)) - 4 * _atan_small(Fraction(1, 239))


def _atan_point(x):
    if x < 0:
        return -_atan_point(-x)
    if x > 1:
        return PI / 2 - _atan_point(1 / x)
    if x <= Fraction(1, 2):
        return _atan_small(x)
    half = Interval(x) / (1 + Interval(1 + x * x).sqrt())
    return 2 * Interval(_atan_small(half.lo).lo, _atan_small(half.hi).hi)


def _acos_interval(cosine):
    if not -1 < cosine.lo <= cosine.hi < 1:
        raise ProofFailure("Angle cosine has not been isolated inside (-1,1)")

    def at(x):
        tangent = ((1 - Interval(x)) / (1 + Interval(x))).sqrt()
        return 2 * Interval(_atan_point(tangent.lo).lo, _atan_point(tangent.hi).hi)

    return Interval(at(cosine.hi).lo, at(cosine.lo).hi)


class ExactField:
    def __init__(self, alpha=None, polynomial=None, interval=None):
        self.alpha = alpha
        self.poly = None
        self.interval = None
        self.root = None
        self._normal = {}
        self._squares = {}
        self._constant_embeddings = {}
        self._domain = None
        if polynomial is not None:
            if alpha is None or interval is None:
                raise ProofFailure("A primitive polynomial needs its symbol and isolating interval")
            self.poly = sp.Poly(polynomial, alpha, domain=sp.QQ)
            if self.poly.degree() < 1 or not self.poly.is_irreducible:
                raise ProofFailure("The primitive polynomial must be irreducible over Q")
            lo, hi = map(sp.Rational, interval)
            if lo == hi and self.poly.degree() == 1 and self.poly.eval(lo) == 0:
                self.interval = Interval(_q(lo))
                self.root = None
                return
            if lo >= hi or self.poly.eval(lo) == 0 or self.poly.eval(hi) == 0:
                raise ProofFailure("The primitive interval must have distinct non-root endpoints")
            if self.poly.count_roots(lo, hi) != 1:
                raise ProofFailure("The primitive interval does not contain exactly one real root")
            self.interval = Interval(_q(lo), _q(hi))
            index = int(self.poly.count_roots(-sp.oo, lo))
            self.root = sp.CRootOf(self.poly.as_expr(), index)

    def normal(self, expr):
        expr = sp.sympify(expr)
        if expr in self._normal:
            return self._normal[expr]
        value = sp.cancel(expr)
        if self.alpha is None or self.alpha not in value.free_symbols:
            value = self._embed_coefficient(value)
        if self.poly is not None and self.alpha in value.free_symbols:
            num, den = sp.fraction(value)
            try:
                def embed_coefficients(polynomial):
                    p = sp.Poly(polynomial, self.alpha, extension=True)
                    return sum(self._embed_coefficient(coefficient) * self.alpha ** power[0]
                               for power, coefficient in p.terms())

                num = embed_coefficients(num)
                den = embed_coefficients(den)
                num = sp.rem(num, self.poly.as_expr(), self.alpha, extension=True)
                den = sp.rem(den, self.poly.as_expr(), self.alpha, extension=True)
                inv = sp.invert(den, self.poly.as_expr(), self.alpha, extension=True)
                value = sp.rem(num * inv, self.poly.as_expr(), self.alpha, extension=True)
            except (sp.PolynomialError, sp.NotInvertible) as exc:
                raise ProofFailure(f"Expression is not a rational function in the primitive field: {expr}") from exc
        value = sp.expand(value)
        self._normal[expr] = value
        return value

    def _embed_coefficient(self, coefficient):
        if coefficient.is_Rational or coefficient.free_symbols:
            return coefficient
        if coefficient in self._constant_embeddings:
            return self._constant_embeddings[coefficient]
        rational, remainder = coefficient.as_coeff_Mul(rational=True)
        if rational != 1:
            value = rational * self._embed_coefficient(remainder)
            self._constant_embeddings[coefficient] = value
            return value
        if coefficient.is_Add:
            value = sum(self._embed_coefficient(term) for term in coefficient.args)
            self._constant_embeddings[coefficient] = value
            return value
        value = sp.simplify(coefficient)
        if self.root is not None and not value.is_Rational:
            try:
                value = sp.to_number_field(value, self.root).as_expr().xreplace({self.root: self.alpha})
            except (sp.IsomorphismFailed, sp.NotAlgebraic, NotImplementedError):
                pass
        self._constant_embeddings[coefficient] = value
        return value

    def zero(self, expr):
        return self.normal(expr) == 0

    def refine(self):
        if self.poly is None:
            return False
        lo, hi = self.interval.lo, self.interval.hi
        if lo == hi:
            return False
        mid = (lo + hi) / 2
        if self.poly.eval(_sq(mid)) == 0:
            # An irreducible rational polynomial with a rational root is linear.
            width = (hi - lo) / 8
            self.interval = Interval(mid - width, mid + width)
        elif self.poly.count_roots(_sq(lo), _sq(mid)) == 1:
            self.interval = Interval(lo, mid)
        else:
            self.interval = Interval(mid, hi)
        return True

    def enclosure(self, expr):
        expr = sp.sympify(expr)
        if expr == self.alpha:
            return self.interval
        if expr.is_Rational:
            return Interval(_q(expr))
        if expr.is_Add:
            answer = Interval(0)
            for arg in expr.args:
                answer += self.enclosure(arg)
            return answer
        if expr.is_Mul:
            answer = Interval(1)
            for arg in expr.args:
                answer *= self.enclosure(arg)
            return answer
        if expr.is_Pow:
            base, exponent = expr.args
            if exponent.is_Integer:
                return self.enclosure(base) ** int(exponent)
            if exponent.is_Rational and exponent.q == 2:
                return self.enclosure(base).sqrt() ** int(exponent.p)
        raise ProofFailure(f"Unsupported expression in rational interval arithmetic: {expr}")

    def sign(self, expr, *, strict=False):
        expr = self.normal(expr)
        if expr == 0:
            if strict:
                raise ProofFailure("A required strictly positive expression is zero")
            return 0
        for _ in range(384):
            try:
                iv = self.enclosure(expr)
                if iv.lo > 0:
                    return 1
                if iv.hi < 0:
                    return -1
            except ProofFailure:
                pass
            if not self.refine():
                break
        raise ProofFailure(f"Could not isolate the sign of {expr}")

    def positive(self, expr):
        if self.sign(expr, strict=True) != 1:
            raise ProofFailure(f"A required positive expression is negative: {expr}")

    def square_root_in_field(self, value):
        """Return a positive base-field root, or None if it is not a square."""
        value = self.normal(value)
        if value in self._squares:
            return self._squares[value]
        if value == 0:
            return sp.S.Zero
        z = sp.Dummy("square_root")
        expression = value
        if self.root is not None:
            expression = expression.subs(self.alpha, self.root)
        extensions = [self.root] if self.root is not None else []
        extensions += [p for p in expression.atoms(sp.Pow)
                       if p.exp.is_Rational and p.exp.q > 1 and not p.free_symbols]
        try:
            value_domain = sp.Poly(value, self.alpha, extension=True).domain if self.root is not None else None
            if self.root is not None and (value_domain.is_QQ or value_domain.is_ZZ):
                if self._domain is None:
                    self._domain = sp.QQ.algebraic_field(self.root)
                coefficients = sp.Poly(value, self.alpha, domain=sp.QQ).all_coeffs()
                element = self._domain.new([sp.QQ.convert(c) for c in coefficients])
                poly = sp.Poly.from_dict({(2,): self._domain.one, (0,): -element}, z,
                                         domain=self._domain)
            elif extensions:
                poly = sp.Poly(z * z - expression, z, extension=extensions)
            else:
                poly = sp.Poly(z * z - expression, z, domain=sp.QQ)
            result = None
            for factor, _ in poly.factor_list()[1]:
                if factor.degree() == 1:
                    candidate = -factor.nth(0) / factor.nth(1)
                    if self.root is not None:
                        candidate = candidate.xreplace({self.root: self.alpha})
                    candidate = self.normal(candidate)
                    if self.sign(candidate) < 0:
                        candidate = -candidate
                    if not self.zero(candidate * candidate - value):
                        raise ProofFailure("Square-root factor did not pass exact substitution")
                    result = candidate
                    break
        except (sp.PolynomialError, sp.CoercionFailed, NotImplementedError) as exc:
            raise ProofFailure(f"Base-field square-root test failed: {exc}") from exc
        self._squares[value] = result
        return result


class RadicalAlgebra:
    """Sparse multiquadratic arithmetic with exact square-class recognition."""

    def __init__(self, field, max_generators=8):
        self.field = field
        self.generators = []
        self.max_generators = max_generators
        self.roots = {}

    def clean(self, terms):
        result = {}
        for mask, value in terms.items():
            value = self.field.normal(value)
            if value != 0:
                result[mask] = value
        return result

    def add(self, a, b):
        result = dict(a)
        for mask, value in b.items():
            result[mask] = result.get(mask, 0) + value
        return self.clean(result)

    def negate(self, a):
        return {mask: -value for mask, value in a.items()}

    def multiply(self, a, b):
        result = defaultdict(lambda: sp.S.Zero)
        for left, x in a.items():
            for right, y in b.items():
                coeff = x * y
                overlap = left & right
                for i, radicand in enumerate(self.generators):
                    if overlap & (1 << i):
                        coeff *= radicand
                result[left ^ right] += coeff
        return self.clean(result)

    def sqrt(self, radicand):
        radicand = self.field.normal(radicand)
        if radicand in self.roots:
            return self.roots[radicand]
        self.field.positive(radicand)
        for mask in range(1 << len(self.generators)):
            product = sp.S.One
            for i, value in enumerate(self.generators):
                if mask & (1 << i):
                    product *= value
            coefficient = self.field.square_root_in_field(radicand / product)
            if coefficient is not None:
                answer = {mask: coefficient}
                self.roots[radicand] = answer
                return answer
        if len(self.generators) >= self.max_generators:
            raise ProofFailure("Fan radical field exceeded the configured generator bound")
        mask = 1 << len(self.generators)
        self.generators.append(radicand)
        answer = {mask: sp.S.One}
        self.roots[radicand] = answer
        return answer


def _topology(faces, vertex_count):
    if vertex_count < 4 or not faces:
        raise ProofFailure("A closed three-dimensional body needs a nonempty sphere complex")
    edges = defaultdict(list)
    links = defaultdict(list)
    vertices = set()
    for face_id, face in enumerate(faces):
        if len(face) not in (3, 4, 5, 6, 8, 10) or len(set(face)) != len(face):
            raise ProofFailure("A face is not a simple permitted regular polygon")
        for i, vertex in enumerate(face):
            vertices.add(vertex)
            previous, following = face[i - 1], face[(i + 1) % len(face)]
            edges[tuple(sorted((vertex, following)))].append(face_id)
            links[vertex].append((previous, following))
    if vertices != set(range(vertex_count)):
        raise ProofFailure("Vertex labels must be exactly 0,...,V-1")
    if any(len(incident) != 2 for incident in edges.values()):
        raise ProofFailure("Every edge must have exactly two incident faces")
    if vertex_count - len(edges) + len(faces) != 2:
        raise ProofFailure("The closed complex does not have Euler characteristic two")
    graph = defaultdict(set)
    for i, j in edges:
        graph[i].add(j)
        graph[j].add(i)

    def connected(adjacency):
        start = next(iter(adjacency))
        seen, pending = {start}, [start]
        while pending:
            for neighbor in adjacency[pending.pop()]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    pending.append(neighbor)
        return len(seen) == len(adjacency)

    if not connected(graph):
        raise ProofFailure("The complex is disconnected")
    for vertex, pairs in links.items():
        link = defaultdict(list)
        for i, j in pairs:
            link[i].append(j)
            link[j].append(i)
        if any(len(adjacent) != 2 for adjacent in link.values()) or not connected(link):
            raise ProofFailure(f"The link of vertex {vertex} is not a circle")
        angle_sum_over_pi = sum(sp.Rational(len(face) - 2, len(face))
                                for face in faces if vertex in face)
        if angle_sum_over_pi >= 2:
            raise ProofFailure(f"Vertex {vertex} does not have strictly positive intrinsic curvature")
    return edges


def _face_coordinates(n):
    # Keeping all y-coordinates as multiples of the same sine avoids artificial
    # square-root extensions introduced by separate trigonometric quotients.
    cosine = sp.cos(2 * sp.pi / n)
    sine = sp.sin(2 * sp.pi / n)
    result = [(sp.S.Zero, sp.S.Zero)]
    x, y = sp.S.Zero, sp.S.Zero
    for step in range(n - 1):
        x += sp.chebyshevt(step, cosine)
        if step:
            y += sine * sp.chebyshevu(step - 1, cosine)
        result.append((sp.simplify(x), sp.simplify(y)))
    return result


def verify_solution(faces, groups, squared_radii, *, primitive_symbol=None,
                    primitive_polynomial=None, isolating_interval=None,
                    primitive_interval=None, equations=(),
                    auxiliary_assignment=None, auxiliary_assignments=None,
                    system=None, max_radical_generators=8):
    """Return CERTIFIED only after all exact geometric obligations pass.

    ``groups`` is a list of vertex lists, and ``squared_radii`` has one exact
    expression per group.  The primitive polynomial must be irreducible over Q.
    Supplied equations/auxiliary assignments are accepted for API compatibility;
    the certificate is independently rebuilt from the face incidences and radii,
    so the solver's particular elimination equations are never trusted.
    """
    report = {"status": "UNVERIFIED", "certified": False, "verified": False,
              "method": "exact radial closure, rational interval winding, linear convex hinges",
              "checks": {}, "limitations": []}
    try:
        if isolating_interval is None:
            isolating_interval = primitive_interval
        field = ExactField(primitive_symbol, primitive_polynomial, isolating_interval)
        groups = [list(map(int, group)) for group in groups]
        faces = [list(map(int, face)) for face in faces]
        if len(groups) != len(squared_radii):
            raise ProofFailure("The number of squared radii differs from the number of groups")
        flat_vertices = [v for group in groups for v in group]
        if len(set(flat_vertices)) != len(flat_vertices):
            raise ProofFailure("Sphere groups overlap")
        vertex_count = len(flat_vertices)
        edges = _topology(faces, vertex_count)
        radii = {}
        for group, radius in zip(groups, squared_radii):
            radius = field.normal(radius)
            field.positive(radius)
            for vertex in group:
                radii[vertex] = radius
        if set(radii) != set(range(vertex_count)):
            raise ProofFailure("Sphere groups do not partition the vertices")
        report["checks"]["combinatorial_sphere"] = True
        report["checks"]["positive_intrinsic_vertex_curvature"] = vertex_count
        report["checks"]["positive_squared_radii"] = True

        edge_d = {}
        for i, j in edges:
            edge_d[i, j] = field.normal(4 * radii[i] * radii[j] - (radii[i] + radii[j] - 1) ** 2)
            field.positive(edge_d[i, j])

        face_projections = {}
        flatness_count = 0
        for face_id, face in enumerate(faces):
            xy = _face_coordinates(len(face))
            for shift in range(len(face)):
                oriented = face[shift:] + face[:shift]
                i, j, k = oriented[:3]
                x, y = xy[2]
                d2 = sp.simplify(x * x + y * y)
                L = field.normal((radii[i] + d2 - (1 + radii[i] - radii[j]) * x - radii[k]) / y)
                for vertex, (xk, yk) in zip(oriented[3:], xy[3:]):
                    residual = radii[i] + xk * xk + yk * yk - (1 + radii[i] - radii[j]) * xk - L * yk - radii[vertex]
                    if not field.zero(residual):
                        raise ProofFailure(f"Face {face_id} fails an exact regular-face flatness equation")
                    flatness_count += 1
                key = tuple(sorted((i, j)))
                field.positive(edge_d[key] - L * L)
                face_projections[face_id, key] = L
        report["checks"]["regular_face_flatness"] = flatness_count
        report["checks"]["strictly_positive_face_cones"] = len(faces)

        coplanar = []
        for key, incident in edges.items():
            hinge = field.normal(sum(face_projections[f, key] for f in incident))
            sign = field.sign(hinge)
            if sign < 0:
                raise ProofFailure(f"Edge {key} has a nonconvex hinge")
            if sign == 0:
                coplanar.append(list(key))
        report["checks"]["convex_hinges"] = len(edges)
        report["coplanar_edges"] = coplanar

        fans = defaultdict(list)
        corner_cache = {}
        for face in faces:
            n = len(face)
            opposite_d2 = sp.simplify(4 * sp.cos(sp.pi / n) ** 2)
            for index, i in enumerate(face):
                j, k = face[index - 1], face[(index + 1) % n]
                a, b, c = radii[i], radii[j], radii[k]
                key = (a, *sorted((b, c), key=sp.default_sort_key), n)
                if key not in corner_cache:
                    Dij = field.normal(4 * a * b - (a + b - 1) ** 2)
                    Dik = field.normal(4 * a * c - (a + c - 1) ** 2)
                    N = field.normal(2 * a * (b + c - opposite_d2) - (a + b - 1) * (a + c - 1))
                    M = field.normal(Dij * Dik)
                    S = field.normal(M - N * N)
                    field.positive(S)
                    corner_cache[key] = (N, M, S)
                fans[i].append(key)

        fan_cache = {}
        fan_reports = []
        for vertex in range(vertex_count):
            keys = fans[vertex]
            signature = tuple(sorted(Counter(keys).items(), key=lambda item: sp.default_sort_key(item[0])))
            if signature in fan_cache:
                fan_reports.append({"vertex": vertex, "same_fan_as": fan_cache[signature]})
                continue
            algebra = RadicalAlgebra(field, max_radical_generators)
            real, imag = {0: sp.S.One}, {}
            target = sp.S.One
            incident_edges = {tuple(sorted((vertex, v))) for face in faces if vertex in face
                              for v in (face[face.index(vertex) - 1], face[(face.index(vertex) + 1) % len(face)])}
            for key in incident_edges:
                target *= edge_d[key]
            target = field.normal(target)
            for key in keys:
                N, M, S = corner_cache[key]
                root = algebra.sqrt(S)
                new_real = algebra.add(algebra.multiply(real, {0: N}),
                                       algebra.negate(algebra.multiply(imag, root)))
                new_imag = algebra.add(algebra.multiply(imag, {0: N}),
                                       algebra.multiply(real, root))
                real, imag = new_real, new_imag
            if algebra.add(real, {0: -target}) or imag:
                raise ProofFailure(f"The exact positive-radical fan at vertex {vertex} does not close")
            winding = None
            for _ in range(256):
                try:
                    total = Interval(0)
                    for key in keys:
                        N, M, S = corner_cache[key]
                        total += _acos_interval(field.enclosure(N) / field.enclosure(M).sqrt())
                    if total.lo > PI.hi and total.hi < 3 * PI.lo:
                        winding = total
                        break
                except ProofFailure:
                    pass
                if not field.refine():
                    break
            if winding is None:
                raise ProofFailure(f"Could not certify one full turn at vertex {vertex}")
            fan_cache[signature] = vertex
            fan_reports.append({"vertex": vertex, "exact_product_identity": True,
                                "angle_sum_interval": winding.strings(),
                                "winding_number": 1,
                                "independent_radicals": len(algebra.generators)})
        report["checks"]["exact_closed_fans"] = vertex_count
        report["checks"]["winding_number_one"] = vertex_count
        report["fan_certificates"] = fan_reports
        report["status"] = "CERTIFIED"
        report["certified"] = True
        report["verified"] = True
        report["conclusion"] = "The proposed squared radii admit an embedded closed convex realization with the specified regular faces and a strictly interior common center."
        if field.interval is not None:
            report["refined_primitive_interval"] = field.interval.strings()
    except (ProofFailure, ValueError, TypeError, NotImplementedError) as exc:
        report["reason"] = str(exc)
    return report
