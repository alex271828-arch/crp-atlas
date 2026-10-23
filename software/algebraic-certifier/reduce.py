"""Bounded exact elimination for the radii encoder.

Numerical values only choose factors and real roots of equations obtained from
the metric.  They are never converted into guessed algebraic constants.  A
successful result is a candidate, to be checked independently against the full
metric and convex realization, not an impossibility or minimality certificate.
"""
from __future__ import annotations

import contextlib
import math
import signal
import threading
import time

import sympy as sp


class ReductionTimeout(Exception):
    pass


@contextlib.contextmanager
def _time_limit(seconds):
    """Interrupt SymPy too, while preserving an enclosing wall-clock timer."""
    if seconds is None or threading.current_thread() is not threading.main_thread():
        yield
        return
    seconds = max(0.001, float(seconds))
    previous_handler = signal.getsignal(signal.SIGALRM)
    previous_timer = signal.getitimer(signal.ITIMER_REAL)
    started = time.monotonic()
    if previous_timer[0] > 0:
        seconds = min(seconds, previous_timer[0])

    def expired(signum, frame):
        raise ReductionTimeout("exact reduction time budget exhausted")

    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        if previous_timer[0] > 0:
            left = max(0.001, previous_timer[0] - (time.monotonic() - started))
            signal.setitimer(signal.ITIMER_REAL, left, previous_timer[1])


def _expr(e):
    if isinstance(e, dict):
        e = e.get("equation", e.get("expression", e.get("polynomial")))
    if isinstance(e, sp.Poly):
        return e.as_expr()
    if isinstance(e, sp.Equality):
        return e.lhs - e.rhs
    return sp.sympify(e)


def _numerator(e):
    e = sp.sympify(e)
    if not e.free_symbols:
        # Constants are "polynomials" to is_polynomial(), but expanding them
        # can leave an exact zero hidden behind algebraic denominators.  A
        # number-field check also handles identities between algebraic powers
        # and isolated roots that rational cancellation alone cannot settle.
        e = sp.cancel(e)
        if (e != 0 and e.is_zero is not False
                and (e.is_algebraic is True or e.has(sp.CRootOf))):
            try:
                if sp.to_number_field(e).coeffs() == [0]:
                    return sp.S.Zero
            except (sp.polys.polyerrors.NotAlgebraic,
                    sp.polys.polyerrors.IsomorphismFailed, NotImplementedError):
                # An unsupported identity stays unresolved; never infer zero
                # from the numerical witness or a floating-point tolerance.
                pass
        return e.as_numer_denom()[0].expand()
    if e.is_polynomial():
        return e.expand()
    return sp.cancel(e).as_numer_denom()[0].expand()


def _value(e, witness):
    try:
        z = complex(e.evalf(25, subs=witness))
        return z.real if abs(z.imag) < 1e-8 else math.nan
    except (TypeError, ValueError, OverflowError):
        return math.nan


def _scaled_residual(e, witness):
    """Coefficient-aware residual, used for branch selection only."""
    variables = sorted(e.free_symbols, key=str)
    val = _value(e, witness)
    if not math.isfinite(val):
        return math.inf
    try:
        p = sp.Poly(e, *variables)
        scale = 0.0
        for powers, coefficient in p.terms():
            term = abs(float(coefficient.evalf()))
            for v, power in zip(variables, powers):
                term *= max(abs(float(witness[v])), 1.0) ** power
            scale += term
        return abs(val) / max(1.0, scale)
    except (KeyError, TypeError, ValueError, sp.PolynomialError):
        return abs(val)


def _select_factors(e, witness, trace, tolerance=1e-6):
    """Select a witnessed factor, retaining ambiguity as a product.

    Removing other branches cannot certify nonexistence: this module only
    proposes solutions, and the independent verifier checks every proposal.
    """
    if e.is_Mul:
        supplied_factors = []
        for term in e.args:
            base, exponent = term.as_base_exp()
            if base.free_symbols and exponent.is_Integer and exponent > 0:
                supplied_factors.append(base)
            elif term.free_symbols:
                supplied_factors.append(term)
        possible = [f for f in supplied_factors
                    if _scaled_residual(f, witness) <= tolerance]
        if possible and len(possible) < len(supplied_factors):
            selected = sp.Mul(*possible)
            trace.append({"operation": "witnessed_factor_selection",
                          "source": e, "selected": selected,
                          "selection_is_not_a_nonexistence_proof": True})
            e = selected
    e = _numerator(e)
    if not e or not e.free_symbols:
        return e
    operation_count = sp.count_ops(e)
    # After rational pivots, cleared denominators can make a bivariate
    # polynomial look large although most of it factors away immediately.
    if operation_count > (20000 if len(e.free_symbols) <= 2 else 1200):
        return e
    try:
        coefficient, factors = sp.factor_list(e, extension=True)
    except (sp.PolynomialError, NotImplementedError, ValueError):
        return e
    if not factors:
        return e
    possible = [f for f, multiplicity in factors
                if _scaled_residual(f, witness) <= tolerance]
    if not possible:
        # A noisy seed must not invent a different equation.
        return sp.Mul(*(f for f, multiplicity in factors)).expand()
    selected = sp.Mul(*possible).expand()
    if len(possible) != len(factors) or any(m != 1 for f, m in factors):
        trace.append({"operation": "witnessed_factor_selection",
                      "source": e, "selected": selected,
                      "selection_is_not_a_nonexistence_proof": True})
    return selected


def _unique(equations):
    answer, seen = [], set()
    for e in equations:
        e = _numerator(e)
        if not e:
            continue
        variables = sorted(e.free_symbols, key=str)
        if variables:
            try:
                e = sp.Poly(e, *variables, extension=True).monic().as_expr()
            except (sp.PolynomialError, ValueError, NotImplementedError):
                pass
        if e not in seen:
            seen.add(e)
            answer.append(e)
    return answer


def _substitute(e, substitutions):
    # Values are kept closed under later substitutions by _eliminate_linear.
    if len(substitutions) == 1:
        v, replacement = next(iter(substitutions.items()))
        numerator, denominator = sp.fraction(replacement)
        if denominator.free_symbols and v in e.free_symbols:
            # Homogeneous Horner substitution clears one common denominator.
            # Substituting a rational expression term by term can instead
            # construct thousands of fractions before cancel() sees them.
            p = sp.Poly(e, v)
            coefficients = p.all_coeffs()
            answer = coefficients[0]
            for i, coefficient in enumerate(coefficients[1:], 1):
                answer = answer * numerator + coefficient * denominator**i
            return answer.expand()
    return _numerator(e.xreplace(substitutions))


def _eliminate_linear(equations, variables, witness, trace, *, protected=()):
    """Cheap linear pivots, including equations linear in one variable."""
    equations = _unique(equations)
    substitutions = {}
    variables = list(variables)
    for iteration in range(len(variables) + 1):
        equations = _unique(_select_factors(e, witness, trace) for e in equations)
        pivots = []
        for index, e in enumerate(equations):
            active = sorted(e.free_symbols.intersection(variables), key=str)
            if not active:
                continue
            if sp.count_ops(e) > 250:
                continue
            for v in active:
                if v in protected:
                    continue
                p = sp.Poly(e, v)
                if p.degree() != 1:
                    continue
                a, b = p.all_coeffs()
                av = _value(a, witness)
                if not math.isfinite(av) or abs(av) < 1e-9:
                    continue
                replacement = sp.cancel(-b / a)
                score = (bool(a.free_symbols), len(replacement.free_symbols),
                         sp.count_ops(replacement), index, str(v))
                pivots.append((score, v, replacement, e, a))
        if not pivots:
            break
        _, v, replacement, source, denominator = min(pivots, key=lambda z: z[0])
        substitutions = {k: sp.cancel(val.subs(v, replacement))
                         for k, val in substitutions.items()}
        substitutions[v] = replacement
        trace.append({"operation": "linear_elimination", "variable": v,
                      "expression": replacement, "source": source,
                      "required_nonzero": denominator})
        equations = [_substitute(e, {v: replacement}) for e in equations]
        variables.remove(v)
    return _unique(equations), variables, substitutions


def _real_roots_near(poly, variable, approximation, tolerance=1e-5):
    """Isolate witnessed roots while preserving the coefficient embedding.

    For coefficients in an algebraic number field, the field norm gives a
    rational polynomial containing the original roots and their coefficient
    conjugates. Exact substitution distinguishes those embeddings; numerical
    proximity selects and orders candidates without proving an identity.
    Rational inputs keep their direct path.
    """
    p = sp.Poly(poly, variable, extension=True)
    rational_coefficients = p.domain.is_QQ or p.domain.is_ZZ
    try:
        rational = p if rational_coefficients else p.norm()
    except (NotImplementedError, sp.PolynomialError):
        return []
    rational = sp.Poly(rational.as_expr(), variable, domain=sp.QQ).sqf_part()
    width = tolerance * max(1.0, abs(approximation))
    left = sp.Rational(str(approximation - width))
    right = sp.Rational(str(approximation + width))
    roots = []
    # Huge coefficient norms can make an entirely wrong factor have a tiny
    # normalized floating residual.  Sturm root counts are exact and do not
    # suffer that problem; isolate only factors meeting the witness interval.
    for factor, multiplicity in rational.factor_list()[1]:
        if factor.count_roots(left, right):
            roots.extend(factor.real_roots(radicals=False))
    scored = []
    for root in roots:
        distance = abs(float(root.evalf(20)) - approximation)
        if distance <= tolerance * max(1.0, abs(approximation)):
            if not rational_coefficients:
                # A tiny coefficient can make a wrong conjugate have an
                # arbitrarily small residual. Its exact algebraic value,
                # rather than a floating tolerance, decides this check.
                residual = p.as_expr().subs(variable, root)
                if residual != 0:
                    try:
                        exact_zero = sp.to_number_field(residual).coeffs() == [0]
                    except (sp.polys.polyerrors.IsomorphismFailed,
                            sp.polys.polyerrors.NotAlgebraic) as exc:
                        # An unsupported exact embedding is unresolved. It
                        # must never be replaced by a numerical zero test.
                        raise ValueError(
                            "Could not verify the algebraic coefficient embedding exactly: "
                            + type(exc).__name__) from exc
                    if not exact_zero:
                        continue
            scored.append((distance, root))
    return [root for _, root in sorted(scored, key=lambda item: item[0])]


def _triangular_solutions(basis, variables, witness, trace):
    """Extract witnessed real points from a triangular zero-dimensional basis."""
    polynomials = [p.as_expr() if isinstance(p, sp.Poly) else p for p in basis]
    if not variables:
        return [{}] if all(_numerator(p) == 0 for p in polynomials) else []
    univariates = [(e, next(iter(e.free_symbols))) for e in polynomials
                  if len(e.free_symbols) == 1]
    if not univariates:
        return []
    polynomial, alpha_var = min(univariates,
                                key=lambda pair: sp.degree(pair[0], pair[1]))
    roots = _real_roots_near(polynomial, alpha_var, witness[alpha_var])
    solutions = []
    for root in roots[:4]:
        expressions = {alpha_var: alpha_var}
        unknowns = set(variables) - {alpha_var}
        for _ in range(len(variables)):
            progress = False
            for e in polynomials:
                e = _substitute(e, {v: ex for v, ex in expressions.items()
                                    if v != alpha_var})
                active = e.free_symbols & unknowns
                if len(active) != 1:
                    continue
                v = next(iter(active))
                p = sp.Poly(e, v)
                if p.degree() != 1:
                    continue
                a, b = p.all_coeffs()
                if abs(_value(a, {alpha_var: float(root.evalf())})) < 1e-12:
                    continue
                expressions[v] = sp.cancel(-b / a)
                unknowns.remove(v)
                progress = True
            if not progress:
                break
        if unknowns:
            # This bounded prototype requires shape-position linear back
            # substitution; it does not mistake a remaining fiber for a root.
            continue
        if isinstance(root, sp.CRootOf):
            modulus = sp.Poly(root.poly.as_expr().subs(root.poly.gen, alpha_var),
                              alpha_var, domain=sp.QQ)
            for v, expression in list(expressions.items()):
                numerator, denominator = sp.cancel(expression).as_numer_denom()
                try:
                    a = sp.Poly(numerator, alpha_var, domain=sp.QQ)
                    b = sp.Poly(denominator, alpha_var, domain=sp.QQ)
                    expressions[v] = (a.rem(modulus) * sp.invert(b, modulus)).rem(modulus).as_expr()
                except (sp.PolynomialError, ValueError, sp.polys.polyerrors.CoercionFailed):
                    pass
        assignment = {v: ex.subs(alpha_var, root) for v, ex in expressions.items()}
        if all(abs(_value(assignment[v], {} ) - witness[v]) <=
               1e-5 * max(1.0, abs(witness[v])) for v in variables):
            trace.append({"operation": "exact_real_root_isolation",
                          "polynomial": polynomial, "variable": alpha_var,
                          "root": root})
            solutions.append(assignment)
    return solutions


def _structural_capacity(equations, variables):
    """Maximum matching in the equation/variable incidence graph.

    A deficient matching prevents a nonempty zero-dimensional complex system:
    some subset of variables occurs in fewer equations than its cardinality.
    This is a cheap reason to defer, never a proof of non-realizability.
    """
    owners = {}
    supports = [e.free_symbols.intersection(variables) for e in equations]

    def augment(variable, visited):
        for i, support in enumerate(supports):
            if i in visited or variable not in support:
                continue
            visited.add(i)
            if i not in owners or augment(owners[i], visited):
                owners[i] = variable
                return True
        return False

    return sum(augment(v, set()) for v in variables)


def _apply_closed(equation, substitutions):
    for variable, replacement in substitutions.items():
        equation = _substitute(equation, {variable: replacement})
    return equation


def _combine_substitutions(first, second):
    return {**{v: sp.cancel(e.xreplace(second)) for v, e in first.items()}, **second}


def _quadratic_pair_consequences(equations, variables, trace, *, limit=24):
    """Cancel matching quadratic leading terms to expose linear pivots.

    These are simple ideal combinations, not resultants: if f=a*v²+... and
    g=b*v²+..., then b*f-a*g is at most linear in v.  Bounds keep this a local
    simplification rather than an uncontrolled polynomial completion.
    """
    small = [e for e in equations if e.free_symbols and sp.count_ops(e) <= 90]
    additions, seen = [], set(equations)
    for variable in variables:
        quadratics = [(e, sp.Poly(e, variable)) for e in small if variable in e.free_symbols]
        quadratics = [(e, p) for e, p in quadratics if p.degree() == 2]
        for i, (first, a) in enumerate(quadratics):
            for second, b in quadratics[i+1:]:
                remainder = sp.expand(b.LC()*first-a.LC()*second)
                if not remainder or remainder in seen or sp.count_ops(remainder) > 120:
                    continue
                if sp.Poly(remainder, variable).degree() > 1:
                    continue
                seen.add(remainder)
                additions.append(remainder)
                trace.append({"operation": "quadratic_leading_term_cancellation",
                              "variable": variable, "first": first, "second": second,
                              "consequence": remainder})
                if len(additions) >= limit:
                    return additions
    return additions


def _staged_linear(equations, variables, witness, trace):
    """Let small patch relations simplify large fan consequences first."""
    equations = list(equations)
    equations += _quadratic_pair_consequences(equations, variables, trace)
    cheap, other = [], []
    for e in equations:
        degree = sp.Poly(e, *sorted(e.free_symbols, key=str)).total_degree() if e.free_symbols else 0
        (cheap if degree <= 3 and sp.count_ops(e) <= 120 else other).append(e)
    small, remaining, first = _eliminate_linear(cheap, variables, witness, trace)
    trace.append({"operation": "low_degree_prepass", "input_equations": len(cheap),
                  "eliminated_variables": tuple(first), "deferred_equations": len(other)})
    other = [_apply_closed(e, first) for e in other]
    combined = [*small, *other]
    combined += _quadratic_pair_consequences(combined, remaining, trace)
    reduced, remaining, second = _eliminate_linear(combined, remaining, witness, trace)
    return reduced, remaining, _combine_substitutions(first, second)


def solve_low_degree(equations, variables, witness, *, max_seconds=None,
                     protected=(), max_variables=7):
    """Close a small exact metric system without body-specific algebra.

    A geometric constructor may supply coordinates, Gram entries, heights,
    or radii as variables.  Linear distance differences and quadratic norm
    cancellations are exhausted before any bounded Groebner step.  Numerical
    values select only candidate factors and embeddings.  The caller must
    independently verify the resulting realization.

    ``solutions`` contains assignments for *all* input variables, including
    those removed by rational pivots.  A partially simplified family is kept
    in the other result fields when no isolated candidate can be extracted.
    """
    started = time.monotonic()
    trace, substitutions = [], {}
    original_variables = tuple(variables)
    remaining = list(variables)
    reduced = list(map(_expr, equations))
    result = {"status": "UNRESOLVED", "solutions": [], "derivation": trace}
    try:
        with _time_limit(max_seconds):
            reduced = _unique(reduced)
            for iteration in range(min(8, len(original_variables) + 1)):
                before = tuple(remaining)
                reduced, remaining, current = _eliminate_linear(
                    reduced, remaining, witness, trace, protected=protected)
                substitutions = _combine_substitutions(substitutions, current)
                trace.append({"operation": "local_metric_linear_pass",
                              "iteration": iteration, "remaining_variables": len(remaining),
                              "remaining_equations": len(reduced)})
                additions = _quadratic_pair_consequences(reduced, remaining, trace)
                if additions:
                    combined = _unique([*reduced, *additions])
                    if len(combined) == len(reduced) and before == tuple(remaining):
                        break
                    reduced = combined
                elif before == tuple(remaining) or not remaining:
                    break
            raw_solutions = _groebner_solutions(
                reduced, remaining, witness, trace, max_variables=max_variables)
            for solution in raw_solutions:
                closed = dict(solution)
                closed.update({v: sp.cancel(e.subs(solution, simultaneous=True))
                               for v, e in substitutions.items()})
                if any(v not in closed or closed[v].free_symbols
                       for v in original_variables):
                    continue
                result["solutions"].append({v: closed[v] for v in original_variables})
            if result["solutions"]:
                result["status"] = "ALGEBRAIC_CANDIDATE"
            else:
                result["reason"] = "no isolated witnessed point in the bounded local metric reduction"
    except ReductionTimeout as exc:
        result.update(reason=str(exc), timed_out=True)
    except (ValueError, TypeError, NotImplementedError, ZeroDivisionError,
            sp.PolynomialError, sp.polys.polyerrors.CoercionFailed,
            sp.polys.polyerrors.NotInvertible) as exc:
        result["reason"] = f"{type(exc).__name__}: {exc}"
    result.update(remaining_equations=reduced, remaining_variables=remaining,
                  substitutions=substitutions, seconds=time.monotonic() - started)
    return result


def _radical_square_classes(model, radius_substitutions, witness, trace):
    """Factor exact radicands and identify equal square classes.

    The returned positive auxiliary is shared only when the radicands differ
    by an exact rational-function square.  Floating data chooses the sign of
    that rational factor; the final verifier independently checks the branch.
    """
    classes, replacements, equations = {}, {}, []
    for corner in model.get("corner_types", ()):
        y = corner["y"]
        radicand = sp.cancel(corner["S"].xreplace(radius_substitutions))
        if sp.count_ops(radicand) > 1500:
            continue
        numerator, denominator = sp.fraction(radicand)
        coefficient, powers = sp.S.One, {}
        for expression, sign in ((numerator, 1), (denominator, -1)):
            symbols = sorted(expression.free_symbols, key=str)
            if not symbols:
                # A fully substituted radicand is an algebraic coefficient,
                # not a polynomial with an inferred generator.
                coefficient *= expression**sign
                continue
            constant, factors = sp.factor_list(expression, *symbols, extension=True)
            coefficient *= constant**sign
            for factor, power in factors:
                symbols = sorted(factor.free_symbols, key=str)
                polynomial = sp.Poly(factor, *symbols, extension=True)
                coefficient *= polynomial.LC()**(power*sign)
                factor = polynomial.monic().as_expr()
                powers[factor] = powers.get(factor, 0) + power*sign
        sign = 1 if _value(coefficient, witness) > 0 else -1
        outside = sp.sqrt(sign*coefficient)
        core = sp.Integer(sign)
        for factor, power in powers.items():
            outside *= factor**(power//2)
            if power % 2:
                core *= factor
        outside, core = sp.cancel(outside), sp.expand(core)
        outside_value, core_value = _value(outside, witness), _value(core, witness)
        if not math.isfinite(outside_value) or not math.isfinite(core_value) or core_value <= 0:
            continue
        if outside_value < 0:
            outside = -outside
        if not core.free_symbols:
            replacements[y] = outside*sp.sqrt(core)
            continue
        if core not in classes:
            auxiliary = sp.Symbol("qr" + str(len(classes)), positive=True)
            classes[core] = auxiliary
            witness[auxiliary] = math.sqrt(core_value)
            equations.append(auxiliary**2-core)
        replacements[y] = outside*classes[core]
    trace.append({"operation": "exact_radical_square_classes",
                  "original_radicals": len(model.get("corner_types", ())),
                  "square_classes": len(classes), "substitution": replacements})
    return replacements, equations, tuple(classes.values())


def _groebner_solutions(equations, variables, witness, trace, *, max_variables=7):
    equations = [_expr(e) for e in equations]
    equations = [e if e.free_symbols else _numerator(e) for e in equations]
    active = set().union(*(e.free_symbols for e in equations)) if equations else set()
    if set(variables) - active:
        trace.append({"operation": "groebner_skipped", "reason": "unconstrained_variables",
                      "variables": tuple(set(variables)-active)})
        return []
    variables = [v for v in variables if v in active]
    if not variables:
        return [{}] if all(e == 0 for e in equations) else []
    if len(variables) > max_variables or len(equations) > 120:
        trace.append({"operation": "groebner_skipped", "reason": "size_bound",
                      "variables": len(variables), "equations": len(equations)})
        return []
    capacity = _structural_capacity(equations, variables)
    if capacity < len(variables):
        trace.append({"operation": "groebner_skipped", "reason": "structurally_underconstrained",
                      "matching": capacity, "variables": len(variables),
                      "equations": len(equations)})
        return []
    if any(not e.free_symbols and e != 0 for e in equations):
        return []
    trace.append({"operation": "groebner_start", "variables": tuple(variables),
                  "equations": len(equations), "input_polynomials": tuple(equations)})
    started = time.monotonic()
    basis = sp.groebner(equations, *variables, order="grevlex", extension=True)
    if not basis.is_zero_dimensional:
        trace.append({"operation": "groebner_incomplete",
                      "reason": "not_zero_dimensional",
                      "seconds": time.monotonic() - started})
        return []
    basis = basis.fglm("lex")
    trace.append({"operation": "groebner_complete", "basis": tuple(basis.exprs),
                  "seconds": time.monotonic() - started})
    return _triangular_solutions(basis.polys, variables, witness, trace)


def _root_field_expressions(radius_values, root, alpha, polynomial):
    """Embed algebraic coefficients into an already isolated root's field.

    Closure equations often leave a geometric parameter together with face
    constants such as sqrt(5).  Those constants may already belong to the
    parameter's field.  Embed their joint coefficient field once, then perform
    rational polynomial arithmetic, rather than choosing a new primitive sum
    of the radii.  to_number_field retains the actual real embedding.
    """
    pairs = [sp.cancel(value.xreplace({root: alpha})).as_numer_denom()
             for value in radius_values]
    polynomials, options = sp.parallel_poly_from_expr(
        [part for pair in pairs for part in pair], alpha, extension=True)
    domain = options["domain"]
    if domain.is_QQ:
        converted = [poly.set_domain(sp.QQ).rem(polynomial)
                     for poly in polynomials]
    elif domain.is_AlgebraicField:
        generator = sp.Poly(
            sp.to_number_field(domain.ext, root).as_expr(alpha),
            alpha, domain=sp.QQ).rem(polynomial)
        zero = sp.Poly(0, alpha, domain=sp.QQ)
        variable = sp.Poly(alpha, alpha, domain=sp.QQ)
        coefficient_cache = {}

        def embed_coefficient(coefficient):
            powers = tuple(coefficient.to_list())
            if powers not in coefficient_cache:
                embedded = zero
                for value in powers:
                    rational = sp.Poly(domain.dom.to_sympy(value),
                                       alpha, domain=sp.QQ)
                    embedded = (embedded * generator + rational).rem(polynomial)
                coefficient_cache[powers] = embedded
            return coefficient_cache[powers]

        converted = []
        for poly in polynomials:
            embedded = zero
            for coefficient in poly.rep.to_list():
                embedded = (embedded * variable
                            + embed_coefficient(coefficient)).rem(polynomial)
            converted.append(embedded)
    else:
        raise sp.polys.polyerrors.CoercionFailed(
            "closure coefficients do not form an algebraic number field")
    return [(converted[i] * sp.invert(converted[i + 1], polynomial))
            .rem(polynomial).as_expr()
            for i in range(0, len(converted), 2)]


def _radii_generate_root_field(expressions, alpha, polynomial):
    """Prove the radii generate the entire root field by rational linear algebra.

    Build the span of 1 and products of radii, reducing modulo the root's
    polynomial.  Closure under multiplication by every radius gives the
    generated Q-algebra, which is a field in this finite extension.  Full rank
    proves that preserving the closure parameter does not enlarge the radius
    field; a closed proper span requires the usual primitive-element path.
    """
    degree = polynomial.degree()
    generators = [sp.Poly(value, alpha, domain=sp.QQ).rem(polynomial)
                  for value in expressions]
    pivots = {}
    basis = []

    def admit(value):
        vector = [value.nth(i) for i in range(degree)]
        for index, row in sorted(pivots.items()):
            scale = vector[index]
            if scale:
                vector = [a - scale*b for a, b in zip(vector, row)]
        for index, coefficient in enumerate(vector):
            if coefficient:
                vector = [entry/coefficient for entry in vector]
                pivots[index] = vector
                basis.append(sp.Poly.from_dict(
                    {(i,): entry for i, entry in enumerate(vector) if entry},
                    alpha, domain=sp.QQ))
                return

    admit(sp.Poly(1, alpha, domain=sp.QQ))
    position = 0
    while position < len(basis):
        if len(basis) == degree:
            return True
        element = basis[position]
        position += 1
        for generator in generators:
            admit((element * generator).rem(polynomial))
            if len(basis) == degree:
                return True
    return False


def _primitive_candidate(radius_values, radius_symbols, model, trace):
    alpha = sp.Symbol("alpha", real=True)
    irrational = [v for v in radius_values if v.is_Rational is not True]
    primitive_value = sp.Integer(0)
    if irrational:
        roots = set().union(*(v.atoms(sp.CRootOf) for v in radius_values))
        expressions = None
        if len(roots) == 1:
            root = next(iter(roots))
            p = sp.Poly(root.poly.as_expr().subs(root.poly.gen, alpha), alpha,
                        domain=sp.QQ)
            try:
                expressions = []
                for value in radius_values:
                    numerator, denominator = sp.cancel(value.xreplace({root: alpha})).as_numer_denom()
                    a = sp.Poly(numerator, alpha, domain=sp.QQ)
                    b = sp.Poly(denominator, alpha, domain=sp.QQ)
                    expressions.append((a.rem(p) * sp.invert(b, p)).rem(p).as_expr())
                polynomial = p.clear_denoms()[1].primitive()[1].as_expr()
                primitive_value = root
            except (sp.PolynomialError, ValueError, sp.polys.polyerrors.CoercionFailed):
                expressions = None
            if expressions is None and threading.current_thread() is threading.main_thread():
                # Optional field simplification must not consume the whole
                # solve budget.  Failure leaves the existing compositum path
                # available; no numerical coefficient recognition is used.
                embedding_started = time.monotonic()
                try:
                    with _time_limit(2.0):
                        expressions = _root_field_expressions(
                            radius_values, root, alpha, p)
                        if not _radii_generate_root_field(expressions, alpha, p):
                            raise ValueError("radii generate a proper subfield of the closure root")
                    polynomial = p.clear_denoms()[1].primitive()[1].as_expr()
                    primitive_value = root
                    trace.append({"operation": "primitive_root_coefficient_embedding",
                                  "degree": p.degree(),
                                  "seconds": time.monotonic() - embedding_started})
                except (ReductionTimeout, sp.polys.polyerrors.BasePolynomialError,
                        ValueError, NotImplementedError):
                    expressions = None
        if expressions is None:
            field = sp.to_number_field(irrational)
            polynomial = field.minpoly.as_expr().subs(field.minpoly.gen, alpha)
            polynomial = sp.Poly(polynomial, alpha).clear_denoms()[1].primitive()[1].as_expr()
            expressions = [sp.to_number_field(v, field).as_expr(alpha)
                           for v in radius_values]
            primitive_value = field.root
        approx = float(primitive_value.evalf(30))
        intervals = sp.Poly(polynomial, alpha).intervals(eps=sp.Rational(1, 10**14))
        matching = [(a, b) for (a, b), multiplicity in intervals
                    if float(a) - 1e-12 <= approx <= float(b) + 1e-12]
        if len(matching) != 1:
            raise ValueError("could not isolate the primitive real embedding")
        interval = matching[0]
    else:
        polynomial = alpha
        expressions = radius_values
        interval = (sp.Integer(0), sp.Integer(0))
    assignments = dict(zip(radius_symbols, expressions))
    auxiliaries = {}
    for corner in model.get("corner_types", ()):
        if not isinstance(corner, dict) or "y" not in corner or "S" not in corner:
            continue
        auxiliaries[corner["y"]] = sp.sqrt(sp.cancel(corner["S"].subs(assignments)))
    return {"status": "ALGEBRAIC_CANDIDATE", "primitive_symbol": alpha,
            "primitive_polynomial": polynomial, "primitive_interval": interval,
            "isolating_interval": interval, "squared_radii": expressions,
            "radius_expressions": assignments,
            "auxiliary_assignments": auxiliaries,
            "auxiliary_assignment": auxiliaries,
            "algebraic_values": dict(zip(radius_symbols, radius_values)),
            "primitive_value": primitive_value,
            "derivation": list(trace)}


def solve_model(model, *, max_seconds=None, cheap_only=False):
    """Return exact radius candidates obtained from the encoder's equations.

    The witness chooses a component only.  No unresolved or rejected component
    is evidence of non-realizability of the body.
    """
    start = time.monotonic()
    trace = []
    result = {"status": "UNRESOLVED", "candidates": [], "derivation": trace,
              "phases": []}
    try:
        with _time_limit(max_seconds):
            radii = tuple(model.get("radius_variables", model.get("radii_symbols", ())))
            if not radii:
                raise ValueError("model has no radius_variables")
            witness = dict(model["witness"])
            # Differences expose the otherwise repeated common radial offset.
            differences = sp.symbols("dr1:" + str(len(radii)), real=True)
            transform = {radii[i]: radii[0] + sum(differences[:i])
                         for i in range(1, len(radii))}
            radius_unknowns = [radii[0], *differences]
            for i, d in enumerate(differences, 1):
                witness[d] = witness[radii[i]] - witness[radii[i - 1]]
            trace.append({"operation": "squared_radius_differences",
                          "substitution": transform})
            phase_start = time.monotonic()
            result["phase"] = "patch_and_radius_relations"
            trace.append({"operation": "phase_start", "phase": result["phase"]})
            flatness = list(model.get("flatness_equations", ()))
            patches = list(model.get("patch_equations", ()))
            consequences = list(model.get("consequence_equations", ()))
            equations = [_substitute(_select_factors(_expr(e), witness, trace), transform)
                         for e in [*flatness, *patches, *consequences]]
            equations = [e for e in equations if e.free_symbols.issubset(radius_unknowns)]
            equations, remaining, eliminated = _staged_linear(
                equations, radius_unknowns, witness, trace)
            solutions = _groebner_solutions(equations, remaining, witness, trace)
            result["phases"].append({"name": result["phase"], "seconds": time.monotonic()-phase_start,
                                     "remaining_variables": len(remaining),
                                     "remaining_equations": len(equations),
                                     "solutions": len(solutions)})
            if not solutions and not cheap_only and model.get("deferred_relations"):
                # A relation is translated only after the cheap metric
                # pivots have supplied a parameterized radius family.
                # Commit completed reductions even if they do not yet isolate
                # a point. An interrupted pass retains the last complete one.
                from encoder import materialize_relation, _primitive_equation
                phase_start = time.monotonic()
                result["phase"] = "deferred_local_angle_closure"
                phase_budget = (min(45.0, max(0.01, (max_seconds-(phase_start-start))*0.85))
                                if max_seconds else 45.0)
                pivot_chain = [(entry["variable"], entry["expression"]) for entry in trace
                               if entry["operation"] == "linear_elimination"]
                trace.append({"operation": "phase_start", "phase": result["phase"]})
                try:
                    with _time_limit(phase_budget):
                        pending = list(enumerate(model["deferred_relations"]))
                        original_templates = {}
                        for closure_pass in range(1, len(remaining) + 2):
                            assignments = {v: sp.cancel(transform.get(v, v).xreplace(eliminated))
                                           for v in radii}
                            deferred = list(equations)
                            next_pending = []
                            for index, relation in pending:
                                if index not in original_templates:
                                    original_templates[index] = materialize_relation(
                                        relation, model["corner_types"])
                                materialized, method, _ = original_templates[index]
                                substituted_first = False
                                if not materialized and len(remaining) <= 2:
                                    # A previously unsupported fan can have a
                                    # small closure after another fan removes
                                    # enough parameters. The encoder enforces
                                    # its pairing and multiplicity bounds.
                                    materialized, method, _ = materialize_relation(
                                        relation, model["corner_types"], assignments)
                                    substituted_first = True
                                if not materialized:
                                    next_pending.append((index, relation))
                                for e in materialized:
                                    # Expose positive Gram factors before the
                                    # generic large-expression factoring bound.
                                    e, denominator = _primitive_equation(e)
                                    if denominator.free_symbols:
                                        trace.append({"operation": "local_angle_denominator",
                                                      "required_nonzero": denominator})
                                    # Preserve chronological committed pivots
                                    # when reducing an original-radius template.
                                    e = _select_factors(e, witness, trace)
                                    if not substituted_first:
                                        e = _substitute(e, transform)
                                        for variable, replacement in pivot_chain:
                                            e = _substitute(e, {variable: replacement})
                                            e = _select_factors(e, witness, trace)
                                    if e.free_symbols.issubset(remaining):
                                        e = _select_factors(e, witness, trace)
                                        if e:
                                            deferred.append(e)
                                trace.append({"operation": "local_angle_template_instantiation",
                                              "closure_pass": closure_pass,
                                              "counts": relation["counts"],
                                              "rhs_pi": relation["rhs_pi"], "method": method,
                                              "substituted_radius_family_first": substituted_first,
                                              "radius_family": assignments})
                            previous_count = len(remaining)
                            deferred, deferred_remaining, second = _staged_linear(
                                deferred, remaining, witness, trace)
                            combined = _combine_substitutions(eliminated, second)
                            equations, remaining, eliminated = deferred, deferred_remaining, combined
                            pivot_chain.extend(second.items())
                            trace.append({"operation": "deferred_radius_reduction",
                                          "closure_pass": closure_pass,
                                          "input_variables": previous_count,
                                          "remaining_variables": len(remaining),
                                          "substitutions": second,
                                          "remaining_equations": len(equations)})
                            solutions = _groebner_solutions(equations, remaining, witness, trace)
                            if solutions or not next_pending or len(remaining) >= previous_count:
                                break
                            pending = next_pending
                except ReductionTimeout:
                    trace.append({"operation": "phase_budget_exhausted", "phase": result["phase"],
                                  "budget_seconds": phase_budget})
                result["phases"].append({"name": result["phase"],
                                         "seconds": time.monotonic()-phase_start,
                                         "solutions": len(solutions)})
            if not solutions and not cheap_only:
                short = list(model.get("short_radical_equations", ()))
                # Small legacy models still get their full equations.  Dense
                # many-radical fans are deliberately deferred to a geometric
                # projection model, rather than expanded speculatively.
                if not short and len(model.get("variables", ())) <= 7:
                    short = list(model.get("equations", ()))
                if short:
                    phase_start = time.monotonic()
                    result["phase"] = "compact_signed_radical_relations"
                    trace.append({"operation": "phase_start", "phase": result["phase"]})
                    phase_budget = min(12.0, max(0.01, (max_seconds-(phase_start-start))/3)) if max_seconds else 12.0
                    try:
                        with _time_limit(phase_budget):
                            radius_assignments = {v: sp.cancel(transform.get(v, v).xreplace(eliminated))
                                                  for v in radii}
                            replacements, square_equations, auxiliaries = _radical_square_classes(
                                model, radius_assignments, witness, trace)
                            compact = []
                            for e in short:
                                e = _apply_closed(_expr(e), radius_assignments)
                                e = _apply_closed(e, replacements)
                                if sp.count_ops(e) <= 1800:
                                    compact.append(e)
                            compact += equations + square_equations
                            active = set().union(*(e.free_symbols for e in compact)) if compact else set()
                            all_variables = sorted(active-set(remaining), key=str) + list(remaining)
                            if len(all_variables) <= 10:
                                compact, compact_remaining, second = _staged_linear(
                                    compact, all_variables, witness, trace)
                                compact_solutions = _groebner_solutions(
                                    compact, compact_remaining, witness, trace)
                                if compact_solutions:
                                    solutions = compact_solutions
                                    eliminated = _combine_substitutions(eliminated, second)
                                    remaining = compact_remaining
                            else:
                                trace.append({"operation": "compact_fan_deferred",
                                              "reason": "too_many_square_classes",
                                              "variables": len(all_variables)})
                    except ReductionTimeout:
                        trace.append({"operation": "phase_budget_exhausted", "phase": result["phase"],
                                      "budget_seconds": phase_budget})
                    result["phases"].append({"name": result["phase"], "seconds": time.monotonic()-phase_start,
                                             "solutions": len(solutions)})
                else:
                    trace.append({"operation": "full_fan_deferred", "reason": "no_compact_signed_system",
                                  "original_variables": len(model.get("variables", ()))})
            for solution in solutions:
                closed = dict(solution)
                for v, e in eliminated.items():
                    closed[v] = sp.cancel(e.subs(solution, simultaneous=True))
                radius_values = [sp.cancel(transform.get(v, v).subs(closed,
                                                                  simultaneous=True))
                                 for v in radii]
                if any(v.free_symbols for v in radius_values):
                    continue
                if any(_value(v, {}) <= 0 for v in radius_values):
                    continue
                candidate = _primitive_candidate(radius_values, radii, model, trace)
                result["candidates"].append(candidate)
            if result["candidates"]:
                result["status"] = "ALGEBRAIC_CANDIDATE"
                result["candidate"] = result["candidates"][0]
            else:
                result["reason"] = "no isolated witnessed radius point extracted within supported elimination limits"
    except ReductionTimeout as exc:
        result["reason"] = str(exc)
        result["timed_out"] = True
    except (ValueError, TypeError, NotImplementedError, ZeroDivisionError,
            sp.PolynomialError, sp.polys.polyerrors.CoercionFailed,
            sp.polys.polyerrors.NotInvertible) as exc:
        result["reason"] = f"{type(exc).__name__}: {exc}"
    if result["candidates"]:
        result["status"] = "ALGEBRAIC_CANDIDATE"
        result["candidate"] = result["candidates"][0]
    if result.get("phase") and not any(p["name"] == result["phase"] for p in result["phases"]):
        result["phases"].append({"name": result["phase"],
                                 "seconds": time.monotonic()-phase_start,
                                 "interrupted": True, "reason": result.get("reason")})
    result["seconds"] = time.monotonic() - start
    return result


reduce_system = solve_model
