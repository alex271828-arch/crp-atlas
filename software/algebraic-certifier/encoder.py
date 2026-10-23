"""Exact, combinatorial radius equations for a closed regular-faced surface.

The numerical witness supplies no coefficients.  It is used only to annotate
the positive radical branch.  All equations are necessary conditions; product
closure alone also permits higher winding and must be checked geometrically.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from functools import reduce
from itertools import combinations, product
from math import gcd, isfinite, sqrt
from typing import Any

import sympy as sp


def squared_short_diagonal(n: int) -> sp.Expr:
    """Squared distance between the two neighbors of a regular n-gon corner."""
    if n < 3:
        raise ValueError("A face must have at least three vertices")
    return sp.simplify(2 + 2 * sp.cos(2 * sp.pi / n))


def _primitive_equation(expr: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    """Clear denominators, retaining the required nonvanishing denominator."""
    numerator, denominator = sp.cancel(expr).as_numer_denom()
    if numerator == 0:
        return sp.S.Zero, denominator
    numerator = sp.factor(numerator)
    if numerator.could_extract_minus_sign():
        numerator = -numerator
    # A numerical nonzero multiplier has no bearing on the equation.
    numerator = numerator.as_coeff_Mul()[1]
    return numerator, denominator


def _complex_power(n: sp.Expr, y: sp.Expr, exponent: int):
    """Real/imaginary parts without asking SymPy to expand a complex product."""
    real = sum((-1) ** j * sp.binomial(exponent, 2 * j)
               * n ** (exponent - 2 * j) * y ** (2 * j)
               for j in range(exponent // 2 + 1))
    imag = sum((-1) ** j * sp.binomial(exponent, 2 * j + 1)
               * n ** (exponent - 2 * j - 1) * y ** (2 * j + 1)
               for j in range((exponent + 1) // 2))
    return real, imag


def _fan_product(counts, corner_types):
    real, imag = sp.S.One, sp.S.Zero
    denominator_factors = Counter()
    for count, corner in zip(counts, corner_types):
        if not count:
            continue
        a, b = _complex_power(corner["N"], corner["y"], count)
        real, imag = real * a - imag * b, real * b + imag * a
        denominator_factors[corner["D_left"]] += count
        denominator_factors[corner["D_right"]] += count
    if any(exponent % 2 for exponent in denominator_factors.values()):
        raise ValueError("Vertex fan has an unpaired radial edge")
    rhs = sp.prod(factor ** (exponent // 2)
                  for factor, exponent in denominator_factors.items())
    return real - rhs, imag


def _normalize_relation(counts, rhs_pi):
    counts = tuple(sp.Rational(x) for x in counts)
    nonzero = [x for x in counts if x]
    if not nonzero:
        return None
    scale = sp.ilcm(*[x.q for x in nonzero]) if len(nonzero) > 1 else nonzero[0].q
    integers = [int(x * scale) for x in counts]
    common = reduce(gcd, (abs(x) for x in integers if x))
    scale = sp.Rational(scale, common)
    integers = tuple(int(x * scale) for x in counts)
    rhs_pi = sp.Rational(rhs_pi) * scale
    if next(x for x in integers if x) < 0:
        integers = tuple(-x for x in integers)
        rhs_pi = -rhs_pi
    return integers, rhs_pi


def _chebyshev_numerator(n, m, degree):
    """M^(degree/2) T_degree(N/sqrt(M)), a polynomial in N and M."""
    x = sp.Dummy("x")
    polynomial = sp.Poly(sp.chebyshevt(degree, x), x)
    return sum(coefficient * n ** power[0] * m ** ((degree - power[0]) // 2)
               for power, coefficient in polynomial.terms())


def _denominator_root(corners, powers):
    """A rational positive square root of prod(M_i**powers_i), if visible.

    Each M is a product of two positive radial edge Grams.  Counting those
    factors avoids introducing unnecessary square-root variables.
    """
    factors = Counter()
    for corner, power in zip(corners, powers):
        factors[corner["D_left"]] += power
        factors[corner["D_right"]] += power
    if any(power % 2 for power in factors.values()):
        return None
    return sp.prod(factor ** (power // 2) for factor, power in factors.items())


def _three_angle_gram(counts, rhs_pi, corners):
    """Cosine Gram determinant for three integer-multiple angles.

    Its sign-sensitive unsquared form is used whenever the denominator
    product is visibly a square.  Otherwise the final square is only a
    necessary consequence; the signed positive-radical equation is retained.
    """
    selected = [(c, abs(k)) for c, k in zip(corners, counts) if k]
    cosine = [_chebyshev_numerator(c["N"], c["M"], k) for c, k in selected]
    norms = [c["M"] ** k for c, k in selected]
    w = sum(cosine[i] ** 2 * norms[(i + 1) % 3] * norms[(i + 2) % 3]
            for i in range(3)) - sp.prod(norms)
    root = _denominator_root([c for c, _ in selected], [k for _, k in selected])
    triple = sp.prod(cosine)
    if root is not None:
        # Divide by the proved positive square root of the denominator
        # product *before* expansion.  For an ordinary three-corner cycle
        # this removes half the degree and avoids a large parasitic factor.
        equation = sum(value ** 2 * sp.cancel(root / norm)
                       for value, norm in zip(cosine, norms)) - root
        equation -= 2 * sp.S.NegativeOne ** rhs_pi * triple
        return equation, "three-angle cosine Gram determinant without squaring"
    return w ** 2 - 4 * triple ** 2 * sp.prod(norms), "three-angle cosine Gram determinant with one final square"


def _four_angle_consequence(counts, rhs_pi, corners):
    selected = [(c, 1 if k > 0 else -1) for c, k in zip(corners, counts) if k]
    a, b, c, d = [entry[0] for entry in selected]
    ratio = _denominator_root([a, b, c, d], [1, 1, -1, -1])
    if ratio is None:
        return None
    # cos(ea*alpha+eb*beta) = (-1)^rhs cos(ec*gamma+ed*delta).
    # The pairwise sine signs are retained in the separate mixed equations.
    u = a["N"] * b["N"] - sp.S.NegativeOne ** rhs_pi * ratio * c["N"] * d["N"]
    v = a["S"] * b["S"]
    w = ratio ** 2 * c["S"] * d["S"]
    return (u ** 2 + v - w) ** 2 - 4 * u ** 2 * v


def _radical_power(corner, exponent):
    """Complex power with y² immediately replaced by its exact radicand."""
    n, y, ss = corner["N"], corner["y"], corner["S"]
    power = abs(exponent)
    real = sum((-1) ** j * sp.binomial(power, 2 * j) * n ** (power - 2 * j) * ss ** j
               for j in range(power // 2 + 1))
    imag = y * sum((-1) ** j * sp.binomial(power, 2 * j + 1)
                   * n ** (power - 2 * j - 1) * ss ** j
                   for j in range((power + 1) // 2))
    return real, imag if exponent >= 0 else -imag


def _radical_product(corners, counts):
    real, imag = sp.S.One, sp.S.Zero
    for corner, count in zip(corners, counts):
        if count:
            a, b = _radical_power(corner, count)
            real, imag = real * a - imag * b, real * b + imag * a
    return real, imag


def _short_radical_relation(counts, rhs_pi, corners):
    support = [i for i, k in enumerate(counts) if k]
    if len(support) > 4 or sum(abs(k) for k in counts) > 6:
        return [], [], "mixed relation exceeds construction bound"
    if len(support) == 1:
        i = support[0]
        angle = sp.pi * rhs_pi / counts[i]
        c = corners[i]
        return [sp.sin(angle) * c["N"] - sp.cos(angle) * c["y"]], [], "single-angle positive-radical relation"
    if not rhs_pi.is_Integer:
        return [], [], "nonintegral turn retained without a mixed expansion"
    positive = [i for i in support if counts[i] > 0]
    negative = [i for i in support if counts[i] < 0]
    if positive and negative:
        left = set(positive)
    else:
        # A balanced two-sided equation has fewer radical monomials than
        # blindly multiplying the entire fan.  Ties use canonical indices.
        choices = [set(support[:1]) | set(extra)
                   for size in range(len(support) - 1)
                   for extra in combinations(support[1:], size)]
        left = min(choices, key=lambda side: (
            abs(sum(abs(counts[i]) * (1 if i in side else -1) for i in support)),
            abs(2 * len(side) - len(support)), tuple(sorted(side))))
    left_counts = [k if i in left else 0 for i, k in enumerate(counts)]
    right_counts = [-k if i not in left else 0 for i, k in enumerate(counts)]
    lr, li = _radical_product(corners, left_counts)
    rr, ri = _radical_product(corners, right_counts)
    root = _denominator_root(corners, [abs(a) - abs(b) for a, b in zip(left_counts, right_counts)])
    sigma = sp.S.NegativeOne ** rhs_pi
    if root is not None:
        return [lr - sigma * root * rr, li - sigma * root * ri], [], "balanced positive-radical closure"
    real, imag = _radical_product(corners, counts)
    # Norms are positive, so the argument condition consists of Im=0 and
    # the prescribed real-part sign.  No artificial denominator radical.
    return [imag], [sigma * real], "positive-radical argument and sign"


def _bounded_short_relations(fan_rows, max_combinations=2048):
    """Small exact combinations of fan rows, with a fixed construction budget."""
    rows = sorted(set(tuple(row["counts"]) for row in fan_rows))
    found = {}
    attempted = 0
    exhausted = False
    def consider(indices, coefficients):
        nonlocal attempted, exhausted
        if attempted >= max_combinations:
            exhausted = True
            return
        attempted += 1
        vector = [sum(coefficient * rows[index][j] for index, coefficient in zip(indices, coefficients))
                  for j in range(len(rows[0]))]
        normalized = _normalize_relation(vector, 2 * sum(coefficients))
        if normalized is None:
            return
        counts, rhs = normalized
        if sum(k != 0 for k in counts) <= 4 and sum(abs(k) for k in counts) <= 6:
            found[normalized] = "bounded signed combination of vertex fans"
    for i, j in combinations(range(len(rows)), 2):
        if not any(a and b for a, b in zip(rows[i], rows[j])):
            continue
        for a, b in ((1, -1), (1, -2), (2, -1)):
            consider((i, j), (a, b))
    for indices in combinations(range(len(rows)), 3):
        # All three rows must belong to one overlap component; otherwise a
        # disjoint fan only makes a candidate less sparse.
        overlap = [(i, j) for i, j in combinations(indices, 2)
                   if any(a and b for a, b in zip(rows[i], rows[j]))]
        if len(overlap) < 2:
            continue
        for signs in product((-1, 1), repeat=2):
            consider(indices, (1,) + signs)
        if exhausted:
            break
    return found, dict(attempted_combinations=attempted, accepted_relations=len(found),
                      combination_budget=max_combinations, budget_exhausted=exhausted)


def _relation_consequences(counts, rhs_pi, corners):
    """Small radius-only consequences, never a claim of branch equivalence."""
    support = [(i, abs(c)) for i, c in enumerate(counts) if c]
    if len(support) == 1:
        index, multiplicity = support[0]
        angle = sp.pi * rhs_pi / counts[index]
        cosine = sp.simplify(sp.cos(angle))
        corner = corners[index]
        if cosine == 0:
            equation = corner["N"]
        elif corner["D_left"] == corner["D_right"]:
            equation = corner["N"] - cosine * corner["D_left"]
        else:
            equation = corner["N"] ** 2 - cosine ** 2 * corner["M"]
        return [equation], "single angle fixed exactly", {"cosine": cosine}
    if not rhs_pi.is_Integer:
        return [], "nonintegral half-turn relation retained for exact checking", {}
    if len(support) == 2:
        (i, p), (j, q) = support
        if max(p, q) > 8:
            return [], "large two-angle relation retained without expansion", {}
        a, b = corners[i], corners[j]
        ca = _chebyshev_numerator(a["N"], a["M"], p)
        cb = _chebyshev_numerator(b["N"], b["M"], q)
        da = (a["M"] ** (p // 2) if p % 2 == 0 else
              a["D_left"] ** p if a["D_left"] == a["D_right"] else None)
        db = (b["M"] ** (q // 2) if q % 2 == 0 else
              b["D_left"] ** q if b["D_left"] == b["D_right"] else None)
        if da is not None and db is not None:
            return [ca * db - sp.S.NegativeOne ** rhs_pi * cb * da], "two-angle Chebyshev relation without squaring", {}
        # cos(p*alpha) = (-1)^rhs cos(q*beta), then square once.
        equation = ca ** 2 * b["M"] ** q - cb ** 2 * a["M"] ** p
        return [equation], "two-angle Chebyshev relation", {}
    radius_support = set().union(*(corners[i]["N"].free_symbols | corners[i]["M"].free_symbols
                                  for i, _ in support))
    # A deferred fan can become inexpensive after metric relations leave a
    # small radius family. With two parameters, admit the higher multiplicity
    # case only when paired radial edge factors avoid the final squaring.
    # Keep both the multiplicity and variable bounds; larger fans stay deferred.
    total_multiplicity = sum(k for _, k in support)
    paired_bivariate = (len(support) == 3 and total_multiplicity <= 6
                        and len(radius_support) <= 2
                        and _denominator_root([corners[i] for i, _ in support],
                                              [k for _, k in support]) is not None)
    if len(support) == 3 and (total_multiplicity == 3 or
                              total_multiplicity == 4 and len(radius_support) <= 2 or
                              total_multiplicity <= 6 and len(radius_support) <= 1 or
                              paired_bivariate):
        equation, method = _three_angle_gram(counts, rhs_pi, corners)
        return [equation], method, {}
    if len(support) == 4 and all(k == 1 for _, k in support) and len(radius_support) <= 2:
        equation = _four_angle_consequence(counts, rhs_pi, corners)
        if equation is not None:
            return [equation], "four-angle paired positive-radical norm", {}
    # 2 alpha +/- 2 beta +/- gamma = integer*pi.  If gamma has equal
    # neighboring radii its cosine has no square-root denominator.  Taking
    # a half angle leaves only sqrt(S_alpha*S_beta), eliminated quadratically.
    if len(support) == 3 and sorted(c for _, c in support) == [1, 2, 2]:
        singleton = next(i for i, c in support if c == 1)
        repeated = [i for i, c in support if c == 2]
        g = corners[singleton]
        if g["D_left"] != g["D_right"]:
            return [], "three-angle relation requires another radical", {}
        a, b = (corners[i] for i in repeated)
        u = a["N"] * b["N"]
        v = a["S"] * b["S"]
        r = (g["D_left"] + sp.S.NegativeOne ** rhs_pi * g["N"]) * a["M"] * b["M"] / (2 * g["D_left"])
        equation = (u ** 2 + v - r) ** 2 - 4 * u ** 2 * v
        return [equation], "three-angle half-angle relation", {"U": u, "V": v, "R": r}
    return [], "relation retained without dense expansion", {}


def _unique_equations(expressions):
    seen, result = set(), []
    for expression in expressions:
        if expression == 0:
            continue
        if expression not in seen:
            seen.add(expression)
            result.append(expression)
    return result


def _clean_mixed_equation(expression, radical_symbols):
    numerator, denominator = sp.fraction(sp.together(expression))
    present = [symbol for symbol in radical_symbols if symbol in numerator.free_symbols]
    if not present:
        return _primitive_equation(expression)
    polynomial = sp.Poly(numerator, *present)
    # Factor the small radius coefficient of each radical monomial, rather
    # than running a large multivariate factorization including every y.
    value = sum(sp.factor(coefficient) * sp.prod(symbol ** power for symbol, power in zip(present, monomial))
                for monomial, coefficient in polynomial.terms())
    return value, denominator


def materialize_relation(relation, corner_types, substitutions=None, *, mixed=False):
    """Translate an exact angle relation after earlier geometric reductions.

    In particular, a half-angle norm need not be expanded in all original
    radii.  Applying an already derived radius family first often turns it
    into a much smaller local closure.  This routine introduces no numerical
    coefficients or additional geometric assumptions.
    """
    substitutions = substitutions or {}
    counts, rhs = relation["counts"], relation["rhs_pi"]
    corners = [dict(corner, **{key: sp.factor(sp.cancel(corner[key].xreplace(substitutions)))
                              for key in ("D_left", "D_right", "N", "M", "S")})
               if counts[index] and substitutions else corner
               for index, corner in enumerate(corner_types)]
    if mixed:
        equations, signs, method = _short_radical_relation(counts, rhs, corners)
        return equations, signs, method
    equations, method, detail = _relation_consequences(counts, rhs, corners)
    return equations, method, detail


def build_system(faces, groups, squared_radii_numeric, *, defer_dense=False) -> dict[str, Any]:
    """Build exact equations; ``groups[j]`` lists vertices on sphere ``s_j``.

    The shape must be a closed, edge-to-edge face complex.  Polygon side
    length is one.  Numerical values are squared radii, one per group.
    The intended solution has the sphere center strictly inside the body.
    """
    faces = tuple(tuple(int(v) for v in face) for face in faces)
    groups = tuple(tuple(int(v) for v in group) for group in groups)
    if not faces or not groups or any(not group for group in groups):
        raise ValueError("Nonempty faces and nonempty sphere classes are required")
    vertices = set(v for face in faces for v in face)
    class_of = {}
    for index, group in enumerate(groups):
        for vertex in group:
            if vertex in class_of:
                raise ValueError(f"Vertex {vertex} occurs in more than one sphere class")
            class_of[vertex] = index
    if set(class_of) != vertices:
        raise ValueError("Sphere classes must partition exactly the face vertices")
    numeric = tuple(float(x) for x in squared_radii_numeric)
    if len(numeric) != len(groups) or not all(isfinite(x) for x in numeric):
        raise ValueError("Supply one finite numerical squared radius per class")
    symbols = tuple(sp.Symbol(f"s{i}", positive=True) for i in range(len(groups)))
    edge_faces = defaultdict(list)
    vertex_corners = defaultdict(list)
    face_flatness = []
    polygon_lengths = {}
    corner_keys = set()
    for face_index, face in enumerate(faces):
        n = len(face)
        if n < 3 or len(set(face)) != n:
            raise ValueError(f"Face {face_index} is not a simple polygon")
        diagonal = squared_short_diagonal(n)
        polygon_lengths[n] = diagonal
        for j, vertex in enumerate(face):
            following = face[(j + 1) % n]
            edge_faces[tuple(sorted((vertex, following)))].append(face_index)
            left, right = sorted((class_of[face[j - 1]], class_of[following]))
            key = (class_of[vertex], left, right, n)
            corner_keys.add(key)
            vertex_corners[vertex].append(key)
        if n > 3:
            # Radius squares restricted to a regular face have only Fourier
            # modes 0,+1,-1, equivalently this cyclic order-three recurrence.
            coefficient = diagonal - 1  # 1+2*cos(2*pi/n)
            row = [symbols[class_of[v]] for v in face]
            for j in range(n):
                face_flatness.append(sp.expand(
                    row[(j + 3) % n] - coefficient * row[(j + 2) % n]
                    + coefficient * row[(j + 1) % n] - row[j]))
    if any(len(incident) != 2 for incident in edge_faces.values()):
        raise ValueError("Every edge must belong to exactly two faces")
    # Exclude disjoint cycles in a vertex link: they do not form one fan.
    for vertex in vertices:
        link = defaultdict(list)
        for face in faces:
            if vertex in face:
                j = face.index(vertex)
                u, v = face[j - 1], face[(j + 1) % len(face)]
                link[u].append(v)
                link[v].append(u)
        if any(len(neighbors) != 2 for neighbors in link.values()):
            raise ValueError(f"Vertex {vertex} does not have a cycle link")
        reached, pending = set(), [next(iter(link))]
        while pending:
            u = pending.pop()
            if u not in reached:
                reached.add(u)
                pending.extend(link[u])
        if len(reached) != len(link):
            raise ValueError(f"Vertex {vertex} has disconnected link")

    flatness = _unique_equations([_primitive_equation(x)[0] for x in face_flatness])
    substitution = {}
    if flatness:
        matrix, rhs = sp.linear_eq_to_matrix(flatness, symbols)
        reduced, pivots = matrix.rref()
        for row_index, pivot in enumerate(pivots):
            substitution[symbols[pivot]] = -sum(
                reduced[row_index, j] * symbols[j] for j in range(len(symbols)) if j not in pivots)
    free_symbols = tuple(s for s in symbols if s not in substitution)
    corners, witness, diagnostics = [], dict(zip(symbols, numeric)), []
    positive = list(symbols)
    corner_index = {}
    for index, key in enumerate(sorted(corner_keys)):
        i, j, k, n = key
        a, b, c = symbols[i], symbols[j], symbols[k]
        dl = sp.expand(4 * a * b - (a + b - 1) ** 2)
        dr = sp.expand(4 * a * c - (a + c - 1) ** 2)
        nn = sp.expand(2 * a * (b + c - polygon_lengths[n]) - (a + b - 1) * (a + c - 1))
        mm = dl * dr
        ss = sp.factor(mm - nn ** 2)
        y = sp.Symbol(f"y{index}", positive=True)
        record = dict(index=index, key=key, center_class=i, adjacent_classes=(j, k),
                      face_size=n, squared_corner_diagonal=polygon_lengths[n],
                      D_left=dl, D_right=dr, N=nn, M=mm, S=ss, y=y)
        corners.append(record)
        corner_index[key] = index
        positive.extend((dl, dr, ss, y))
        value = float(ss.evalf(subs=witness))
        if value > 0:
            witness[y] = sqrt(value)
        else:
            diagnostics.append(f"Corner {index} has nonpositive numerical Gram radicand {value}")
    fan_counts_by_vertex = {}
    equivalent_fans = defaultdict(list)
    for vertex in sorted(vertices):
        counts = [0] * len(corners)
        for key in vertex_corners[vertex]:
            counts[corner_index[key]] += 1
        counts = tuple(counts)
        fan_counts_by_vertex[vertex] = counts
        equivalent_fans[counts].append(vertex)
    fan_rows = [dict(vertices=tuple(vs), counts=counts) for counts, vs in equivalent_fans.items()]
    fan_equations = []
    for row in fan_rows:
        row["equations"] = _fan_product(row["counts"], corners)
        fan_equations.extend(row["equations"])

    relations = {}
    def add_relation(counts, rhs, origin):
        normalized = _normalize_relation(counts, rhs)
        if normalized is not None and normalized not in relations:
            relations[normalized] = origin
    for row in fan_rows:
        add_relation(row["counts"], 2, "vertex fan")
    for i, a in enumerate(fan_rows):
        for b in fan_rows[i + 1:]:
            # Disjoint-support differences cannot cancel any angles.
            if any(x and y for x, y in zip(a["counts"], b["counts"])):
                add_relation([x - y for x, y in zip(a["counts"], b["counts"])], 0, "difference of vertex fans")
    # Exact row operations discover relations needing more than one pairwise
    # subtraction.  Keep the augmented constant column to retain the turn.
    augmented = sp.Matrix([list(row["counts"]) + [2] for row in fan_rows])
    rr, _ = augmented.rref()
    inconsistent = False
    for row in rr.tolist():
        if not any(row[:-1]) and row[-1] != 0:
            inconsistent = True
        elif any(row[:-1]):
            add_relation(row[:-1], row[-1], "rational row reduction of vertex fans")
    extra_relations, search_diagnostics = _bounded_short_relations(fan_rows)
    for (counts, rhs), origin in extra_relations.items():
        add_relation(counts, rhs, origin)
    sparse_relations, consequences, denominators = [], [], []
    deferred_relations = []
    short_radicals, short_signs = [], []
    radical_symbols = [corner["y"] for corner in corners]
    method_counts = Counter()
    for (counts, rhs), origin in sorted(relations.items(), key=lambda item: (sum(x != 0 for x in item[0][0]), sum(abs(x) for x in item[0][0]))):
        if defer_dense and (sum(abs(x) for x in counts) >= 5 or
                            sum(x != 0 for x in counts) >= 4):
            record = dict(counts=counts, rhs_pi=rhs, origin=origin,
                          method="deferred local angle closure", equations=[], detail={},
                          short_radical_equations=[], positive_real_parts=[],
                          mixed_method="deferred until geometric substitutions")
            sparse_relations.append(record)
            deferred_relations.append(record)
            method_counts[record["method"]] += 1
            continue
        equations, method, detail = _relation_consequences(counts, rhs, corners)
        cleaned = []
        for equation in equations:
            equation, denominator = _primitive_equation(equation.subs(substitution, simultaneous=True))
            if denominator.free_symbols:
                denominators.append(denominator)
            if equation != 0:
                cleaned.append(equation)
                consequences.append(equation)
        mixed, signs, mixed_method = _short_radical_relation(counts, rhs, corners)
        mixed_cleaned = []
        for equation in mixed:
            equation, denominator = _clean_mixed_equation(
                equation.subs(substitution, simultaneous=True), radical_symbols)
            if denominator.free_symbols:
                denominators.append(denominator)
            if equation != 0:
                mixed_cleaned.append(equation)
                short_radicals.append(equation)
        signs = [sign.subs(substitution, simultaneous=True) for sign in signs]
        short_signs.extend(signs)
        method_counts[method] += 1
        sparse_relations.append(dict(counts=counts, rhs_pi=rhs, origin=origin,
                                     method=method, equations=cleaned, detail=detail,
                                     short_radical_equations=mixed_cleaned,
                                     positive_real_parts=signs, mixed_method=mixed_method))
    if inconsistent:
        consequences.insert(0, sp.S.One)
    quadratics = [c["y"] ** 2 - c["S"] for c in corners]
    equations = flatness + quadratics + fan_equations
    # Original full expressions are already polynomial.  Only the flatness
    # substitution can introduce algebraic denominators; record all of them.
    reduced_equations = []
    for equation in equations:
        expression = equation.subs(substitution, simultaneous=True)
        numerator, denominator = sp.fraction(sp.together(expression))
        if denominator.free_symbols:
            denominators.append(denominator)
        if numerator != 0:
            reduced_equations.append(numerator)
    metadata = dict(coefficient_origin="exact regular-face geometry only",
                    requires_interior_center=True, winding_requires_verification=True,
                    consequence_equations_require_branch_verification=True,
                    denominator_nonzero=_unique_equations(denominators),
                    angle_system_inconsistent=inconsistent, diagnostics=diagnostics,
                    short_angle_search=search_diagnostics,
                    short_angle_methods=dict(method_counts),
                    short_radical_real_parts_positive=_unique_equations(short_signs))
    return dict(faces=faces, groups=groups, numeric_radii=numeric,
                radii_symbols=symbols, radius_variables=symbols,
                variables=symbols + tuple(c["y"] for c in corners),
                class_of=class_of, corner_types=corners,
                corner_index=corner_index, polygon_squared_corner_diagonals=polygon_lengths,
                flatness_equations=flatness, flatness_substitution=substitution,
                free_radii_symbols=free_symbols, fan_rows=fan_rows,
                fan_counts_by_vertex=fan_counts_by_vertex,
                fan_equations=fan_equations, radical_equations=quadratics,
                equations=equations, reduced_equations=_unique_equations(reduced_equations),
                consequence_equations=_unique_equations(consequences),
                short_radical_equations=_unique_equations(short_radicals),
                short_radical_sign_conditions=_unique_equations(short_signs),
                sparse_relations=sparse_relations, deferred_relations=deferred_relations,
                positivity=_unique_equations(positive),
                witness=witness, metadata=metadata)


def serialize_system(system):
    """JSON-compatible inspection record (not an executable expression format)."""
    def convert(value):
        if isinstance(value, sp.Basic):
            return str(value)
        if isinstance(value, dict):
            return {str(key): convert(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [convert(item) for item in value]
        return value
    return convert(system)
