"""Optional, bounded, exactly checked presentations of certified squared radii.

This module changes how an existing certificate is written, never its proof.
It retains rational functions and factors, shares repeated expressions, and
tries elementary alternative generators. No approximate value identifies an
algebraic constant, and no polynomial/Galois-group database is used.
"""
from __future__ import annotations

from fractions import Fraction
import itertools
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from check import expression
from reduce import ReductionTimeout, _time_limit
from verify import Interval


def _field_value(expr, symbol, modulus):
    """Evaluate a rational expression in QQ[symbol]/modulus, exactly."""
    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator = sp.Poly(numerator, symbol, domain=sp.QQ).rem(modulus)
    denominator = sp.Poly(denominator, symbol, domain=sp.QQ).rem(modulus)
    if denominator.is_zero:
        raise ValueError('A presentation denominator vanishes in the source field')
    return (numerator * denominator.invert(modulus)).rem(modulus)


def _fraction_expression(numerator, denominator):
    """Remove common rational scale without expanding the quotient."""
    coefficients = numerator.all_coeffs() + denominator.all_coeffs()
    scale = sp.ilcm(*[coefficient.q for coefficient in coefficients])
    numerator = numerator.mul_ground(scale)
    denominator = denominator.mul_ground(scale)
    common = sp.igcd(*[abs(int(c)) for c in
                       numerator.all_coeffs() + denominator.all_coeffs()])
    if denominator.LC() < 0:
        common = -common
    numerator = numerator.exquo_ground(common)
    denominator = denominator.exquo_ground(common)
    return numerator.as_expr() / denominator.as_expr()


def _rational_forms(poly, modulus):
    """Euclidean rational reconstructions, all exactly congruent to poly.

    Along the extended Euclidean sequence, remainder = multiplier*poly mod f.
    Thus remainder/multiplier is an equivalent short rational expression.
    Unlike choosing one fixed degree bound, considering the whole sequence
    allows the visible formula size to choose its own numerator/denominator
    balance.
    """
    yield poly.as_expr()
    if poly.is_zero:
        return
    old_r, r = modulus, poly
    old_v = sp.Poly(0, *poly.gens, domain=sp.QQ)
    v = sp.Poly(1, *poly.gens, domain=sp.QQ)
    while not r.is_zero:
        if not v.is_zero and sp.gcd(v, modulus).degree() == 0:
            candidate = _fraction_expression(r, v)
            # Retain an explicit exact check, independent of the invariant.
            num, den = sp.fraction(candidate)
            if (sp.Poly(den, *poly.gens, domain=sp.QQ)*poly
                    - sp.Poly(num, *poly.gens, domain=sp.QQ)).rem(modulus).is_zero:
                yield candidate
        quotient, remainder = old_r.div(r)
        old_r, r = r, remainder
        old_v, v = v, old_v - quotient*v


def _short_form(poly, modulus, hints=()):
    """Use only an exactly equivalent shorter expression."""
    best = poly.as_expr()
    best_size = len(str(best))
    candidates = list(hints) + list(_rational_forms(poly, modulus))
    # Cheap textual ranking prevents factorization of enormous losing forms.
    candidates.sort(key=lambda value: len(str(value)))
    for value in candidates:
        if len(str(value)) > max(80, best_size*2):
            continue
        try:
            if _field_value(value-poly.as_expr(), poly.gen, modulus).is_zero:
                numerator, denominator = sp.fraction(value)
                factored_numerator = sp.factor(numerator)
                factored_denominator = sp.factor(denominator)
                for styled in (value, factored_numerator/denominator,
                               numerator/factored_denominator,
                               factored_numerator/factored_denominator):
                    size = len(str(styled))
                    if size < best_size:
                        best, best_size = styled, size
        except (ValueError, sp.PolynomialError, sp.NotInvertible):
            continue
    return best


def _score(polynomial, definitions, radii):
    """Characters in the visible polynomial, shared definitions and radii.

    The exact source map belongs to the audit data and is not printed as a
    radius formula. Q is displayed using the already defined squared radii,
    so expanded Q coefficients cannot inflate the display score.
    """
    return len(str(polynomial)) + sum(len(str(v)) for v in radii) + sum(
        len(str(symbol)) + 3 + len(str(value)) for symbol, value in definitions)


def _share(polynomial, radii):
    definitions, reduced = sp.cse(radii, symbols=sp.numbered_symbols('q'),
                                  order='canonical')
    return _prune_definitions(polynomial, definitions, reduced)


def _prune_definitions(polynomial, definitions, reduced):
    # Raw CSE can invent names for powers used only twice, making a human
    # presentation longer. Inline any such name whose total cost is greater
    # than the text it saves, including its definition and later references.
    changed = True
    while changed:
        changed = False
        old_score = _score(polynomial, definitions, reduced)
        for index in range(len(definitions)-1, -1, -1):
            symbol, value = definitions[index]
            substitution = {symbol: value}
            trial_definitions = definitions[:index] + [
                (name, expr.xreplace(substitution))
                for name, expr in definitions[index+1:]]
            trial_radii = [expr.xreplace(substitution) for expr in reduced]
            if _score(polynomial, trial_definitions, trial_radii) <= old_score:
                definitions, reduced = trial_definitions, trial_radii
                changed = True
                break
    return definitions, reduced


def _affine_coefficients(value, basis):
    """Express a field polynomial in a small existing rational span, if possible."""
    degree = max([value.degree(), *(part.degree() for part in basis)])
    matrix = sp.Matrix([[part.nth(row) for part in basis]
                        for row in range(int(degree)+1)])
    vector = sp.Matrix([value.nth(row) for row in range(int(degree)+1)])
    try:
        answer, parameters = matrix.gauss_jordan_solve(vector)
    except ValueError:
        return None
    # The caller appends a basis element only when it increases the span.
    if parameters.rows:
        return None
    return list(answer)


def _anchored_forms(polynomial, radii, anchor=0):
    """Write radii as a common offset plus affine functions of the generator.

    Divide each radius difference by t in the polynomial power basis.  Short
    rational reconstruction of the quotient exposes such forms as a+b*t+c.
    Exact linear algebra shares quotient coefficients already in the rational
    span of earlier ones.  This is independent of the body's name or field.
    """
    t = polynomial.gen
    variable = sp.Poly(t, t, domain=sp.QQ)
    base = radii[anchor]
    a = sp.Symbol('a')
    definitions = [(a, _short_form(base, polynomial))]
    basis = [sp.Poly(1, t, domain=sp.QQ)]
    names = [sp.S.One]
    forms = []
    for radius in radii:
        quotient, remainder = (radius-base).div(variable)
        coefficients = _affine_coefficients(quotient, basis)
        if coefficients is None:
            name = sp.Symbol('b'+str(len(names)-1))
            definitions.append((name, _short_form(quotient, polynomial)))
            basis.append(quotient)
            names.append(name)
            coefficient = name
        else:
            coefficient = sum(c*name for c, name in zip(coefficients, names))
        forms.append(a + t*coefficient + remainder.as_expr())
    return _prune_definitions(polynomial.as_expr(), definitions, forms)


def _interval_eval(poly, lo, hi):
    value = Interval(0)
    point = Interval(Fraction(lo), Fraction(hi))
    for coefficient in poly.all_coeffs():
        value = value*point + Fraction(coefficient)
    return sp.Rational(value.lo), sp.Rational(value.hi)


def _generator_interval(source_poly, source_interval, generator, target_poly):
    """Enclose the image of the known source root, then count target roots.

    Rational interval evaluation and Sturm root counts choose the embedding;
    no approximate nearest-root match is used.
    """
    lo, hi = source_interval
    if target_poly.degree() == 1:
        root = -target_poly.nth(0)/target_poly.nth(1)
        return [root, root]
    for _ in range(40):
        left, right = _interval_eval(generator, lo, hi)
        if (left < right and target_poly.eval(left) != 0
                and target_poly.eval(right) != 0
                and target_poly.count_roots(left, right) == 1):
            # The interval-evaluation fractions can be hundreds of digits.
            # A coarser rational enclosure is easier to print and has the
            # same rigorous Sturm certificate when it still isolates one root.
            scale = sp.Integer(10)**10
            simple_left = sp.floor(left*scale)/scale
            simple_right = sp.ceiling(right*scale)/scale
            if (simple_left < simple_right and target_poly.eval(simple_left) != 0
                    and target_poly.eval(simple_right) != 0
                    and target_poly.count_roots(simple_left, simple_right) == 1):
                return [simple_left, simple_right]
            return [left, right]
        if lo == hi:
            raise ValueError('A non-linear generator has a rational source root')
        lo, hi = source_poly.refine_root(lo, hi, eps=(hi-lo)/1024)
    raise ValueError('Could not isolate the transformed real generator')


def _change_generator(generator, source_poly, radii, target_symbol):
    """Change a full primitive basis by exact rational linear algebra."""
    n = source_poly.degree()
    powers = [sp.Poly(1, source_poly.gen, domain=sp.QQ)]
    for _ in range(n):
        powers.append((powers[-1]*generator).rem(source_poly))
    basis = DomainMatrix.from_Matrix(sp.Matrix([
        [power.nth(row) for power in powers[:-1]] for row in range(n)]
    )).convert_to(sp.QQ)
    rhs = DomainMatrix.from_Matrix(sp.Matrix([
        [power.nth(row) for power in [powers[-1]] + list(radii)]
        for row in range(n)]
    )).convert_to(sp.QQ)
    try:
        solved = basis.inv().matmul(rhs).to_Matrix()
    except (ValueError, ZeroDivisionError):
        return None
    coefficients = [sum(solved[row, column]*target_symbol**row for row in range(n))
                    for column in range(solved.cols)]
    minimal = sp.Poly(target_symbol**n-coefficients[0], target_symbol, domain=sp.QQ)
    minimal = minimal.clear_denoms()[1].primitive()[1]
    if minimal.LC() < 0:
        minimal = -minimal
    return minimal, [sp.Poly(value, target_symbol, domain=sp.QQ)
                     for value in coefficients[1:]]


def _make_result(poly, interval, source_definition, description,
                 definitions, radii, baseline_score, method, started):
    return {'format_version': 1,
            'generator': {'symbol': str(poly.gen), 'polynomial': str(poly.as_expr()),
                          'isolating_interval': [str(v) for v in interval],
                          'definition_in_source': str(source_definition),
                          'description': description},
            'definitions': [{'symbol': str(symbol), 'expression': str(value)}
                            for symbol, value in definitions],
            'squared_radii': [str(value) for value in radii],
            'Q_factored': '*'.join('(T-('+str(value)+'))' for value in radii),
            'score': _score(poly.as_expr(), definitions, radii),
            'baseline_score': baseline_score, 'method': method, 'verified': True,
            'seconds': time.monotonic()-started}


def _evaluate_field_expression(value, mapping, modulus):
    """Evaluate the expression tree with reduction after every operation.

    Substituting a degree-15 source map into a degree-16 expression and only
    then expanding would create huge intermediate coefficients. Polynomial
    Horner-style field arithmetic avoids that expansion altogether.
    """
    zero = sp.Poly(0, modulus.gen, domain=sp.QQ)
    one = sp.Poly(1, modulus.gen, domain=sp.QQ)
    cache = {}

    def visit(expr):
        if expr in cache:
            return cache[expr]
        if expr.is_Rational:
            result = sp.Poly(expr, modulus.gen, domain=sp.QQ)
        elif expr in mapping:
            result = mapping[expr]
        elif expr.is_Add:
            result = zero
            for term in expr.args:
                result = result + visit(term)
        elif expr.is_Mul:
            result = one
            for term in expr.args:
                result = (result*visit(term)).rem(modulus)
        elif expr.is_Pow and expr.exp.is_Integer:
            base, exponent = visit(expr.base), int(expr.exp)
            if exponent < 0:
                if base.is_zero:
                    raise ValueError('A presentation denominator vanishes')
                base, exponent = base.invert(modulus), -exponent
            result = one
            while exponent:
                if exponent % 2:
                    result = (result*base).rem(modulus)
                exponent //= 2
                if exponent:
                    base = (base*base).rem(modulus)
        else:
            raise ValueError('Presentation is not rational arithmetic in its generator')
        cache[expr] = result.rem(modulus)
        return cache[expr]
    return visit(value)


def _check_candidate(source_poly, source_radii, generator, target_poly,
                     definitions, formulas):
    """Check every presentation expression back in the certified source field."""
    known = {target_poly.gen: generator}
    if not _evaluate_field_expression(target_poly.as_expr(), known, source_poly).is_zero:
        raise ValueError('The generator map does not satisfy its polynomial')
    for symbol, value in definitions:
        known[symbol] = _evaluate_field_expression(value, known, source_poly)
    for old, value in zip(source_radii, formulas):
        mapped = _evaluate_field_expression(value, known, source_poly)
        if mapped != old:
            raise ValueError('A shortened squared radius differs from its source')
    if len(formulas) != len(source_radii):
        raise ValueError('The presentation changed the number of squared radii')


def build_presentation(saved_radii, *, derivation=None, max_seconds=5.0,
                       on_improvement=None):
    """Return the best completed exactly equivalent presentation within budget.

    ``saved_radii`` is the existing radii.json object, with strings parsed by
    the certificate checker's restricted arithmetic parser. The baseline is
    emitted first. Timeouts and unsuitable alternate generators preserve the
    last completed result; they do not affect the original proof or solver
    status. An external process timeout can use ``on_improvement`` to persist
    these same checkpoints.
    """
    started = time.monotonic()
    source_symbol = sp.Symbol(str(saved_radii['primitive_symbol']))
    target_symbol = sp.Symbol('t')
    source_names = {str(source_symbol): source_symbol}
    source_poly = sp.Poly(expression(saved_radii['primitive_polynomial'], source_names),
                          source_symbol, domain=sp.QQ)
    if source_poly.degree() < 1:
        raise ValueError('A positive-degree primitive polynomial is required')
    source_interval = [expression(value, {}) for value in saved_radii['isolating_interval']]
    original = [expression(value, source_names) for value in saved_radii['squared_radii']]
    baseline_poly = sp.Poly(source_poly.as_expr().xreplace({source_symbol: target_symbol}),
                           target_symbol, domain=sp.QQ)
    baseline_radii = [value.xreplace({source_symbol: target_symbol}) for value in original]
    baseline_score = _score(baseline_poly.as_expr(), [], baseline_radii)
    best = _make_result(baseline_poly, source_interval, source_symbol,
                        't = '+str(source_symbol), [], baseline_radii,
                        baseline_score, 'original certificate presentation', started)
    if on_improvement is not None:
        on_improvement(best)
    if max_seconds <= 0:
        return best

    from sympy.polys.matrices.exceptions import DMNonInvertibleMatrixError
    recoverable = (ValueError, sp.PolynomialError, sp.NotInvertible,
                   DMNonInvertibleMatrixError, ZeroDivisionError)

    def accept(poly, candidate_radii, generator, description, hints=()):
        forms = [_short_form(value, poly, hints[i] if i < len(hints) else ())
                 for i, value in enumerate(candidate_radii)]
        definitions, forms = _share(poly.as_expr(), forms)
        publish(poly, generator, description, definitions, forms)
        if len(candidate_radii) > 1:
            for anchor in range(len(candidate_radii)):
                definitions, forms = _anchored_forms(poly, candidate_radii, anchor)
                publish(poly, generator, description, definitions, forms)

    def publish(poly, generator, description, definitions, forms):
        nonlocal best
        score = _score(poly.as_expr(), definitions, forms)
        if score >= best['score']:
            return
        _check_candidate(source_poly, source_radii, generator, poly, definitions, forms)
        interval = (source_interval if generator == identity else
                    _generator_interval(source_poly, source_interval, generator, poly))
        candidate = _make_result(poly, interval, generator.as_expr(), description,
                                 definitions, forms, baseline_score,
                                 'exact rational reconstruction, factorization and shared expressions', started)
        best = candidate
        if on_improvement is not None:
            on_improvement(candidate)

    try:
        remaining = max_seconds-(time.monotonic()-started)
        if remaining <= 0:
            return best
        with _time_limit(remaining):
            source_radii = [_field_value(value, source_symbol, source_poly) for value in original]
            identity = sp.Poly(source_symbol, source_symbol, domain=sp.QQ).rem(source_poly)
            initial = [sp.Poly(value.as_expr().xreplace({source_symbol: target_symbol}),
                               target_symbol, domain=sp.QQ) for value in source_radii]
            # First preserve the original formulas, then reconstruct/factor in
            # the original generator. A degree-one source needs no new basis.
            hints = [[value] for value in baseline_radii]
            for hint_object in (saved_radii.get('presentation_hints'), derivation):
                if not isinstance(hint_object, dict):
                    continue
                raw_hints = hint_object.get('squared_radii', [])
                if isinstance(raw_hints, list) and len(raw_hints) == len(hints):
                    for slot, raw in zip(hints, raw_hints):
                        try:
                            slot.append(expression(raw, source_names).xreplace(
                                {source_symbol: target_symbol}))
                        except (ValueError, SyntaxError, TypeError):
                            pass
            accept(baseline_poly, initial, identity, 't = '+str(source_symbol), hints)
            if source_poly.degree() <= 1:
                return best

            candidates = []
            seen = {identity.as_expr()}
            # Adjacent differences first, then the remaining pairs; no body
            # label or known coordinate formula influences these candidates.
            pairs = list(zip(range(len(source_radii)-1), range(1, len(source_radii))))
            pairs += [pair for pair in itertools.combinations(range(len(source_radii)), 2)
                      if pair not in pairs]
            possibilities = [(source_radii[j]-source_radii[i], f't = s{j} - s{i}')
                             for i, j in pairs]
            possibilities += [(value, f't = s{i}') for i, value in enumerate(source_radii)]
            for generator, description in possibilities:
                if generator.degree() <= 0 or generator.as_expr() in seen:
                    continue
                seen.add(generator.as_expr())
                try:
                    converted = _change_generator(generator, source_poly, source_radii,
                                                   target_symbol)
                    if converted is not None:
                        poly, values = converted
                        # Complete the inexpensive basis changes first. The
                        # shorter polynomials get the presentation budget first.
                        candidates.append((len(str(poly.as_expr())), poly, values,
                                           generator, description))
                except recoverable:
                    continue
            candidates.sort(key=lambda candidate: candidate[0])
            for _, poly, values, generator, description in candidates:
                try:
                    accept(poly, values, generator, description)
                except recoverable:
                    continue
    except ReductionTimeout:
        best = dict(best, budget_exhausted=True)
    except recoverable as exc:
        best = dict(best, stopped_reason=type(exc).__name__)
    best['seconds'] = time.monotonic()-started
    return best
