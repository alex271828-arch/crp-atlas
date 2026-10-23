"""Exact positive-root arithmetic; zero tests never use a numerical tolerance."""

import ast
from functools import lru_cache
import sympy as sp
from vendor.radial import (
    ExactField as InheritedField,
    RadicalAlgebra,
    ProofFailure,
    Interval,
    _topology,
)


def require(condition, message):
    if not condition:
        raise ProofFailure(message)


class ExactField(InheritedField):
    """Treat number-field embedding as a proposal; certify every positive root."""

    def _embed_coefficient(self, coefficient):
        if coefficient.is_Rational or coefficient.free_symbols:
            return coefficient
        if coefficient in self._constant_embeddings:
            return self._constant_embeddings[coefficient]
        if coefficient.is_Add or coefficient.is_Mul:
            value = coefficient.func(*(self._embed_coefficient(x) for x in coefficient.args))
        elif coefficient.is_Pow and coefficient.exp.is_Integer:
            value = self._embed_coefficient(coefficient.base) ** coefficient.exp
        elif coefficient.is_Pow and coefficient.exp.is_Rational and coefficient.exp.q == 2:
            value = coefficient
            if self.root is not None:
                try:
                    proposal = sp.to_number_field(coefficient, self.root).as_expr()
                except (sp.IsomorphismFailed, sp.NotAlgebraic, NotImplementedError):
                    proposal = None
                if proposal is not None:
                    candidate = proposal.xreplace({self.root: self.alpha})
                    target = self._embed_coefficient(coefficient.base) ** int(coefficient.exp.p)
                    require(
                        self.zero(candidate * candidate - target),
                        "Embedding proposal fails exact square identity",
                    )
                    require(
                        self.sign(candidate) == 1,
                        "Embedding proposal selects the wrong real square root",
                    )
                    value = candidate
        else:
            raise ProofFailure("Unsupported exact algebraic coefficient")
        self._constant_embeddings[coefficient] = value
        return value


def parse(text, symbol=None, field=None):
    """Restricted arithmetic, checking divisors before symbolic cancellation."""
    require(type(text) in (str, int), "Exact arithmetic text or integer required")

    def visit(n):
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return sp.Integer(n.value)
        if isinstance(n, ast.Name) and symbol is not None and n.id == str(symbol):
            return symbol
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd, ast.USub)):
            a = visit(n.operand)
            return a if isinstance(n.op, ast.UAdd) else -a
        if isinstance(n, ast.BinOp):
            a, b = visit(n.left), visit(n.right)
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                require(not (field.zero(b) if field else b == 0), "Zero denominator")
                return a / b
            if isinstance(n.op, ast.Pow):
                require(b.is_Integer and abs(b) <= 10000, "Only bounded integer powers allowed")
                if b < 0:
                    require(not (field.zero(a) if field else a == 0), "Zero negative power")
                return a**b
        if (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Name)
            and n.func.id == "sqrt"
            and len(n.args) == 1
            and not n.keywords
        ):
            a = visit(n.args[0])
            require(
                field.sign(a) >= 0 if field else a.is_nonnegative is True,
                "A square root must have a proved nonnegative real radicand",
            )
            return sp.sqrt(a)
        raise ProofFailure("Unsupported arithmetic syntax")

    return visit(ast.parse(str(text), mode="eval").body)


class Algebra:
    def __init__(self, saved):
        self.symbol = sp.Symbol(saved["primitive_symbol"])
        poly = parse(saved["primitive_polynomial"], self.symbol)
        interval = [parse(x) for x in saved["isolating_interval"]]
        self.field = ExactField(self.symbol, poly, interval)
        self.rad = RadicalAlgebra(self.field, 12)
        self.zero = {}
        self.one = {0: sp.S.One}

    def base(self, x):
        x = self.field.normal(x)
        return {0: x} if x != 0 else {}

    @lru_cache(maxsize=65536)
    def read_base(self, x):
        return self.field.normal(parse(x, self.symbol, self.field))

    def add(self, a, b):
        return self.rad.add(a, b)

    def neg(self, a):
        return self.rad.negate(a)

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        return self.rad.multiply(a, b)

    def scale(self, a, b):
        return self.mul(a, self.base(b))

    def sqrt(self, x):
        return self.rad.sqrt(x)

    def total(self, values):
        out = {}
        for x in values:
            out = self.add(out, x)
        return out

    def dot(self, a, b):
        return self.total(self.mul(x, y) for x, y in zip(a, b))

    def vsub(self, a, b):
        return [self.sub(x, y) for x, y in zip(a, b)]

    def cross(self, a, b):
        return [
            self.sub(self.mul(a[i], b[j]), self.mul(a[j], b[i]))
            for i, j in [(1, 2), (2, 0), (0, 1)]
        ]

    def matvec(self, m, v):
        return [self.dot(row, v) for row in m]

    def norm(self, p, g):
        return self.dot(p, self.matvec(g, p))

    def sign(self, a):
        return self._sign(tuple(sorted(self.rad.clean(a).items())))

    @lru_cache(maxsize=65536)
    def _sign(self, items):
        a = dict(items)
        if not a:
            return 0
        for _ in range(384):
            try:
                total = Interval(0)
                for mask, c in a.items():
                    val = self.field.enclosure(c)
                    for i, r in enumerate(self.rad.generators):
                        if mask & (1 << i):
                            val *= self.field.enclosure(r).sqrt()
                    total += val
                if total.lo > 0:
                    return 1
                if total.hi < 0:
                    return -1
            except ProofFailure:
                pass
            if not self.field.refine():
                break
        raise ProofFailure("Could not certify a radical expression sign")

    def encode(self, a):
        result = {}
        for mask, c in sorted(a.items()):
            candidates = [c, sp.factor(c)]
            short = min(candidates, key=lambda x: len(str(x)))
            require(self.field.zero(short - c), "Simplification changed coefficient")
            result[str(mask)] = str(short)
        return result

    def decode(self, data):
        require(isinstance(data, dict), "Sparse radical element required")
        out = {}
        for key, value in data.items():
            require(key.isdigit() and str(int(key)) == key, "Invalid radical mask")
            mask = int(key)
            require(mask < (1 << len(self.rad.generators)), "Radical mask exceeds basis")
            out[mask] = self.read_base(value)
        return self.rad.clean(out)


@lru_cache(maxsize=None)
def chord(n, gap):
    return sp.simplify((sp.sin(sp.pi * gap / n) / sp.sin(sp.pi / n)) ** 2)


def topology(case, saved):
    faces = case["faces"]
    groups = saved["groups"]
    flat = [v for g in groups for v in g]
    require(all(type(v) is int for f in faces for v in f), "Noninteger vertex label")
    require(all(type(v) is int for v in flat), "Noninteger sphere label")
    require(
        len(set(flat)) == len(flat) and set(flat) == set(range(len(flat))),
        "Invalid sphere partition",
    )
    require(groups == case["groups"], "Input and output sphere groups differ")
    _topology(faces, len(flat))
    require(len(groups) == len(saved["squared_radii"]), "Radius count differs")
    return faces, len(flat)
