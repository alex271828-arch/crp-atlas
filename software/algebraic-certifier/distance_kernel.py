"""Exact local Euclidean constructions used by the v3 solver.

The primitive is an intersection of an affine space and one metric sphere.
Numerical witnesses may order pivot attempts; they never supply coefficients,
discard an exact equation, prove a rank, or select a solution branch.

Returned conditions are *obligations*, not certificates.  In particular a
symbolic pivot describes the open branch where that pivot is nonzero.  A caller
must retain this fact and verify its final candidate independently.
"""
from itertools import combinations
import math

import sympy as sp


class ConstructionLimit(ValueError):
    """A configured algebraic size bound prevented a local construction."""


def _normal(expression, max_ops):
    expression = sp.sympify(expression)
    if expression.has(sp.Float):
        raise ValueError("local construction coefficients must be exact; floats belong only in the witness")
    if sp.count_ops(expression) > max_ops:
        raise ConstructionLimit("local expression exceeds the operation bound")
    result = sp.cancel(expression)
    if sp.count_ops(result) > max_ops:
        raise ConstructionLimit("normalized local expression exceeds the operation bound")
    return result


def _vector(value):
    result = sp.Matrix(value)
    if result.cols != 1:
        raise ValueError("a point must be a column vector")
    return result


def _metric(metric, dimension):
    result = sp.eye(dimension) if metric is None else sp.Matrix(metric)
    if result.shape != (dimension, dimension) or result != result.T:
        raise ValueError("metric must be a symmetric square matrix of the point dimension")
    return result


def _unique_conditions(expressions, max_ops):
    result = []
    for expression in expressions:
        expression = _normal(expression, max_ops)
        numerator = expression.as_numer_denom()[0]
        if numerator.is_zero:
            raise ValueError("construction requires an identically zero denominator to be nonzero")
        if numerator.is_number:
            continue
        if numerator.could_extract_minus_sign():
            numerator = -numerator
        if numerator not in result:
            result.append(numerator)
    return result


def _denominators(expressions):
    return [sp.denom(expression) for expression in expressions]


def _numeric_weight(expression, witness):
    if not witness and expression.free_symbols:
        return 0.0
    try:
        value = complex(expression.evalf(16, subs=witness or {}))
        if math.isfinite(value.real) and math.isfinite(value.imag):
            return abs(value)
    except (TypeError, ValueError, OverflowError):
        pass
    return 0.0


def _independent_minor(matrix, witness, max_minors, max_ops):
    """Return a nonzero maximum-size minor; all zero decisions are exact.

    If the bound is exhausted, raise rather than falsely reporting a lower
    rank.  A nonconstant determinant comes with its nonzero branch condition.
    """
    count = 0
    for rank in range(min(matrix.shape), 0, -1):
        candidates = []
        for rows in combinations(range(matrix.rows), rank):
            for columns in combinations(range(matrix.cols), rank):
                count += 1
                if count > max_minors:
                    raise ConstructionLimit("exact local rank search exceeded its minor bound")
                minor = matrix.extract(rows, columns)
                determinant = _normal(minor.det(method="domain-ge"), max_ops)
                if determinant != 0:
                    candidates.append((rows, columns, determinant))
        if candidates:
            # Prefer a pivot nonzero at the guide and then a small expression.
            # No numerical tolerance is used to establish algebraic rank.
            candidates.sort(key=lambda item: (
                _numeric_weight(item[2], witness) < 1e-12,
                sp.count_ops(item[2]),
                -_numeric_weight(item[2], witness), item[0], item[1]))
            return (*candidates[0], count)
    return (), (), sp.S.One, count


def distance_intersection(anchors, squared_distances, *, metric=None,
                          coordinates=None, witness=None, parameter_prefix="t",
                          max_dimension=4, max_minors=256, max_ops=20000):
    """Intersect prescribed metric spheres without eliminating free geometry.

    ``anchors[j]`` is an exact coordinate vector and ``squared_distances[j]``
    is its prescribed squared distance to the unknown point.  The result's
    ``affine_point`` contains ``free_parameters``.  It solves the chosen linear
    equations identically; *all* other affine residuals and ``residual_norm``
    must still vanish.  This handles overdetermined or parameter-dependent data
    without silently asserting that they are consistent.

    ``coordinates`` may name the original unknown coordinates for readable
    equations.  The free parameters are real Dummy symbols, avoiding collisions
    with shape parameters occurring in anchor expressions.
    """
    anchors = [_vector(anchor) for anchor in anchors]
    squared_distances = list(map(sp.sympify, squared_distances))
    if not anchors or len(anchors) != len(squared_distances):
        raise ValueError("supply at least one anchor and one squared distance per anchor")
    dimension = anchors[0].rows
    if not 1 <= dimension <= max_dimension:
        raise ConstructionLimit("point dimension exceeds the local construction bound")
    if any(anchor.rows != dimension for anchor in anchors):
        raise ValueError("all anchors must have the same dimension")
    H = _metric(metric, dimension)
    coordinates = tuple(coordinates) if coordinates is not None else tuple(
        sp.Dummy(f"x{i}", real=True) for i in range(dimension))
    if len(coordinates) != dimension or len(set(coordinates)) != dimension:
        raise ValueError("coordinates must name distinct unknowns of the point dimension")
    if any(not isinstance(value, sp.Symbol) for value in coordinates):
        raise ValueError("coordinates must be symbols")
    x = sp.Matrix(coordinates)
    norm = lambda p: (p.T*H*p)[0]
    a0, ell0 = anchors[0], squared_distances[0]
    rows = [2*(anchor-a0).T*H for anchor in anchors[1:]]
    A = sp.Matrix.vstack(*rows) if rows else sp.zeros(0, dimension)
    b = sp.Matrix([_normal(norm(anchor)-norm(a0)+ell0-ell, max_ops)
                   for anchor, ell in zip(anchors[1:], squared_distances[1:])])
    if not rows:
        b = sp.zeros(0, 1)
    pivot_rows, pivot_columns, determinant, minor_count = _independent_minor(
        A, witness, max_minors, max_ops)
    rank = len(pivot_rows)
    free_columns = tuple(i for i in range(dimension) if i not in pivot_columns)
    parameters = tuple(sp.Dummy(f"{parameter_prefix}{i}", real=True)
                       for i in range(len(free_columns)))
    point = sp.zeros(dimension, 1)
    for column, parameter in zip(free_columns, parameters):
        point[column] = parameter
    if rank:
        block = A.extract(pivot_rows, pivot_columns)
        right = b.extract(pivot_rows, [0])
        if free_columns:
            right -= A.extract(pivot_rows, free_columns)*sp.Matrix(parameters)
        solved = block.inv(method="DM")*right
        for column, value in zip(pivot_columns, solved):
            point[column] = _normal(value, max_ops)
    residuals = [_normal(value, max_ops) for value in A*point-b]
    residuals = [value for value in residuals if value != 0]
    residual_norm = _normal(norm(point-a0)-ell0, max_ops)
    original_equations = [_normal(norm(x-anchor)-ell, max_ops)
                          for anchor, ell in zip(anchors, squared_distances)]
    conditions = _unique_conditions(
        [determinant]+_denominators(list(H)+[v for a in anchors for v in a]
                                    +squared_distances+list(point)+residuals+[residual_norm]),
        max_ops)
    positive_metric = [_normal(H[:i, :i].det(), max_ops) for i in range(1, dimension+1)]
    return {
        "affine_point": point, "free_parameters": parameters,
        "free_columns": free_columns, "coordinates": coordinates,
        "residual_norm": residual_norm, "affine_residuals": residuals,
        "rank": rank, "metric": H, "anchors": anchors,
        "squared_distances": squared_distances,
        "linear_matrix": A, "linear_rhs": b,
        "original_equations": original_equations,
        "nonzero_conditions": conditions,
        "positive_metric_conditions": positive_metric,
        "provenance": {
            "operation": "subtract_distance_equations_then_one_norm",
            "anchor_count": len(anchors), "dimension": dimension,
            "pivot_rows": list(pivot_rows), "pivot_columns": list(pivot_columns),
            "exact_minors_examined": minor_count,
            "rank_is_conditional_on_nonzero_pivot": bool(determinant.free_symbols),
            "witness_use": "pivot_order_only" if witness else "none",
        },
    }


def quadratic_branches(intersection, *, known_point=None, max_ops=20000):
    """Return both exact branches of a line/sphere intersection.

    A known point allows the second solution to be expressed rationally by
    Vieta.  Its validity equations are returned, never assumed.  Branches are
    unfiltered: a witness must not turn one selected candidate into a claimed
    complete set of realizations.
    """
    parameters = intersection["free_parameters"]
    if len(parameters) != 1:
        raise ValueError("quadratic branches require a one-dimensional affine intersection")
    parameter = parameters[0]
    polynomial = sp.Poly(intersection["residual_norm"], parameter)
    if polynomial.degree() != 2:
        raise ValueError("the residual norm is not quadratic in its free parameter")
    a, b, c = polynomial.all_coeffs()
    discriminant = _normal(b*b-4*a*c, max_ops)
    known_conditions = []
    if known_point is None:
        roots = [(-b-sp.sqrt(discriminant))/(2*a),
                 (-b+sp.sqrt(discriminant))/(2*a)]
        method = "quadratic_formula_both_branches"
    else:
        known_point = _vector(known_point)
        if known_point.shape != intersection["affine_point"].shape:
            raise ValueError("known point dimension does not match the intersection")
        root = known_point[intersection["free_columns"][0]]
        known_conditions = [_normal(v, max_ops) for v in
                            intersection["affine_point"].subs(parameter, root)-known_point]
        known_conditions += [_normal(polynomial.as_expr().subs(parameter, root), max_ops)]
        known_conditions += [_normal(e.subs(parameter, root), max_ops)
                             for e in intersection["affine_residuals"]]
        known_conditions = [value for value in known_conditions if value != 0]
        roots = [root, _normal(-b/a-root, max_ops)]
        method = "known_solution_and_vieta_reflection"
    branches = []
    for root in roots:
        point = intersection["affine_point"].subs(parameter, root).applyfunc(
            lambda value: _normal(value, max_ops))
        branches.append({
            "parameter_value": root, "point": point,
            "residual_equations": [_normal(e.subs(parameter, root), max_ops)
                                   for e in intersection["affine_residuals"]],
        })
    return {
        "branches": branches, "discriminant": discriminant,
        "nonnegative_conditions": [discriminant],
        "nonzero_conditions": _unique_conditions(
            intersection["nonzero_conditions"]+[a]
            +_denominators(roots+known_conditions
                           +[value for branch in branches for value in branch["point"]]), max_ops),
        "known_point_equations": known_conditions,
        "positive_metric_conditions": intersection["positive_metric_conditions"],
        "provenance": {"operation": method, "both_branches_retained": True},
    }


def reflect_across_anchor_hyperplane(point, anchors, *, metric=None,
                                      witness=None, max_ops=20000):
    """Reflect a point in the affine span of d independent plane anchors.

    This is the rational second branch of their equal-distance intersection;
    it works in any supported dimension and needs no named polyhedral piece.
    """
    point = _vector(point)
    H = _metric(metric, point.rows)
    anchors = [_vector(anchor) for anchor in anchors]
    distances = [((point-anchor).T*H*(point-anchor))[0] for anchor in anchors]
    intersection = distance_intersection(anchors, distances, metric=H,
                                         witness=witness, max_ops=max_ops)
    if intersection["rank"] != point.rows-1:
        raise ValueError("anchors must span an affine hyperplane")
    result = quadratic_branches(intersection, known_point=point, max_ops=max_ops)
    result["reflected_point"] = result["branches"][1]["point"]
    result["provenance"]["operation"] = "equal_distance_profile_reflection"
    return result


def gram_extension(anchor_gram, inner_products, squared_norm, *,
                   ambient_dimension=3, max_ops=20000):
    """Project a new vector onto independent anchors via a Schur complement.

    ``schur_residual`` is the squared length of its orthogonal component.  It
    must be nonnegative, or zero when the anchors already span the ambient
    space.  Sylvester conditions are returned to make the positive metric
    hypothesis explicit.
    """
    G = sp.Matrix(anchor_gram)
    if G.rows != G.cols or G != G.T or G.rows > ambient_dimension:
        raise ValueError("anchor Gram matrix must be symmetric and of admissible dimension")
    inner = _vector(inner_products) if G.rows else sp.zeros(0, 1)
    if inner.rows != G.rows:
        raise ValueError("one inner product is needed per anchor")
    determinant = _normal(G.det(), max_ops)
    if determinant == 0:
        raise ValueError("Gram anchors must be independent")
    weights = (G.inv(method="DM")*inner).applyfunc(lambda value: _normal(value, max_ops))
    residual = _normal(sp.sympify(squared_norm)-(inner.T*weights)[0], max_ops)
    return {
        "projection_weights": weights, "schur_residual": residual,
        "equations": [residual] if G.rows == ambient_dimension else [],
        "nonnegative_conditions": [] if G.rows == ambient_dimension else [residual],
        "positive_anchor_conditions": [_normal(G[:i, :i].det(), max_ops)
                                       for i in range(1, G.rows+1)],
        "nonzero_conditions": _unique_conditions(
            [determinant]+_denominators(list(G)+list(inner)+list(weights)+[residual]), max_ops),
        "provenance": {"operation": "gram_schur_complement",
                       "anchor_rank": G.rows, "ambient_dimension": ambient_dimension},
    }


def complete_gram_entry(anchor_gram, left_inner, right_inner, left_norm, right_norm,
                        *, unknown=None, ambient_dimension=3, max_ops=20000):
    """Complete one Gram entry with at most one unspanned direction.

    With d-1 anchors, the two residual components lie on a line.  Their dot
    product has the two signs of sqrt(r_left*r_right).  With d anchors it is
    already determined.  All rank and positivity conditions remain explicit.
    """
    G = sp.Matrix(anchor_gram)
    if G.rows not in (ambient_dimension-1, ambient_dimension):
        raise ValueError("Gram completion needs d-1 or d independent anchors")
    left = gram_extension(G, left_inner, left_norm,
                          ambient_dimension=ambient_dimension, max_ops=max_ops)
    right = gram_extension(G, right_inner, right_norm,
                           ambient_dimension=ambient_dimension, max_ops=max_ops)
    projected_dot = _normal((_vector(left_inner).T*right["projection_weights"])[0], max_ops)
    unknown = sp.Dummy("g", real=True) if unknown is None else sp.sympify(unknown)
    if G.rows == ambient_dimension:
        equation = unknown-projected_dot
        branches = [projected_dot]
    else:
        product = _normal(left["schur_residual"]*right["schur_residual"], max_ops)
        equation = _normal((unknown-projected_dot)**2-product, max_ops)
        branches = [projected_dot-sp.sqrt(product), projected_dot+sp.sqrt(product)]
    return {
        "unknown": unknown, "equation": equation, "branches": branches,
        "other_equations": left["equations"]+right["equations"],
        "nonnegative_conditions": left["nonnegative_conditions"]+right["nonnegative_conditions"],
        "positive_anchor_conditions": left["positive_anchor_conditions"],
        "nonzero_conditions": _unique_conditions(
            left["nonzero_conditions"]+right["nonzero_conditions"], max_ops),
        "provenance": {"operation": "rank_constrained_gram_entry",
                       "both_signs_retained": G.rows < ambient_dimension},
    }


def gram_rank_relations(matrix, rank, *, witness=None, max_minors=1024,
                        max_ops=20000):
    """Give fraction-free local equations for a specified PSD Gram rank.

    Choose an exact nonzero principal ``rank`` minor A, then return every
    equation ``det(A)*B_ij - b_i.T*adj(A)*b_j = 0`` in its Schur complement.
    On the explicitly retained branch det(A) != 0 these equations impose
    exactly the requested rank.  Positive principal minors of A then imply
    positive semidefiniteness of the full matrix.  Rank is supplied by the
    caller; it is never inferred as a theorem from approximate coordinates.

    ``witness`` only prioritizes anchor choices.  Its entries may be a mapping
    from symbolic variables to numerical guide values, or a numerical matrix
    of the same shape.  No guide value becomes an equation coefficient.
    """
    matrix = sp.Matrix(matrix)
    if matrix.rows != matrix.cols or matrix != matrix.T:
        raise ValueError("a Gram matrix must be square and symmetric")
    if not isinstance(rank, int) or rank < 0 or rank > matrix.rows:
        raise ValueError("specified Gram rank must be between zero and the matrix size")
    if any(value.has(sp.Float) for value in matrix):
        raise ValueError("Gram matrix coefficients must be exact")
    if rank == 0:
        constraints = [_normal(matrix[i, j], max_ops)
                       for i in range(matrix.rows) for j in range(i, matrix.cols)]
        constraints = list(dict.fromkeys(e for e in constraints if e != 0))
        return {"constraints": constraints, "equations": constraints,
                "anchor_indices": (), "rank": rank,
                "nonzero_conditions": _unique_conditions(_denominators(list(matrix)), max_ops),
                "positive_conditions": [], "psd_conditions": [],
                "provenance": {"operation": "zero_rank_gram_block", "rank": 0}}
    numerical_matrix = None
    numerical_mapping = witness
    if witness is not None and not hasattr(witness, "items"):
        numerical_matrix = sp.Matrix(witness)
        if numerical_matrix.shape != matrix.shape:
            raise ValueError("numerical Gram guide must have the same shape")
        numerical_mapping = None
    candidates = []
    for number, indices in enumerate(combinations(range(matrix.rows), rank), 1):
        if number > max_minors:
            raise ConstructionLimit("Gram anchor search exceeded its principal-minor bound")
        anchor = matrix.extract(indices, indices)
        determinant = _normal(anchor.det(method="domain-ge"), max_ops)
        if determinant == 0:
            continue
        guide = (numerical_matrix.extract(indices, indices).det()
                 if numerical_matrix is not None else determinant)
        weight = _numeric_weight(guide, numerical_mapping)
        candidates.append((weight < 1e-12, sp.count_ops(determinant), -weight,
                           indices, determinant, anchor))
    if not candidates:
        raise ValueError("no nonzero principal anchor minor exists for the requested Gram rank")
    _, _, _, indices, determinant, anchor = min(candidates, key=lambda item: item[:4])
    other = tuple(i for i in range(matrix.rows) if i not in indices)
    adjugate = anchor.adjugate(method="berkowitz")
    constraints = []
    for offset, i in enumerate(other):
        column_i = matrix.extract(indices, [i])
        for j in other[offset:]:
            column_j = matrix.extract(indices, [j])
            relation = _normal(determinant*matrix[i, j]-(column_i.T*adjugate*column_j)[0], max_ops)
            if relation != 0 and relation not in constraints:
                constraints.append(relation)
    psd = [_normal(anchor[:i, :i].det(), max_ops) for i in range(1, rank+1)]
    conditions = _unique_conditions([determinant]+_denominators(list(matrix)+constraints), max_ops)
    return {"constraints": constraints, "equations": constraints,
            "anchor_indices": indices, "rank": rank,
            "nonzero_conditions": conditions,
            "positive_conditions": psd, "psd_conditions": psd,
            "anchor_determinant": determinant,
            "provenance": {"operation": "fraction_free_gram_schur_complement",
                           "rank": rank, "anchor_indices": list(indices),
                           "witness_use": "anchor_order_only" if witness is not None else "none",
                           "conditions_are_obligations": True}}
