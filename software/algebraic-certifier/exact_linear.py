"""Linear and polynomial arithmetic over explicit algebraic coefficient fields.

Generic expression-domain arithmetic can keep algebraically equal coefficients
in different forms and cause large intermediate expressions.  These helpers
collect coefficients into one exact number field before doing arithmetic.  The
polynomial variables are always explicit; they never become field constants.
"""
from __future__ import annotations

import sympy as sp
from sympy.matrices.exceptions import NonInvertibleMatrixError, NonSquareMatrixError
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.exceptions import DMNonInvertibleMatrixError


def _exact_expressions(expressions, variables):
    allowed = set(variables)
    values = [sp.sympify(expression) for expression in expressions]
    for value in values:
        if not isinstance(value, sp.Expr):
            raise TypeError("Expected scalar polynomial expressions")
        if value.has(sp.Float):
            raise ValueError("Exact coefficient arithmetic does not accept floats")
        if value.free_symbols - allowed:
            raise ValueError("All polynomial variables must be supplied explicitly")
    return values


def _algebraic_field(domain):
    field = domain.get_field()
    if not (field.is_QQ or field.is_AlgebraicField or field.is_GaussianField):
        raise ValueError("Coefficients must belong to an exact algebraic number field")
    return field


def exact_inverse(matrix):
    """Invert a constant algebraic matrix, returning a regular SymPy Matrix.

    Rational and algebraic coefficients are supported, with no numerical
    recognition or expression-domain fallback.  Non-square and singular inputs
    raise the same exception classes as ``Matrix.inv``.
    """
    matrix = sp.Matrix(matrix)
    if matrix.rows != matrix.cols:
        raise NonSquareMatrixError("A matrix must be square to invert")
    _exact_expressions(matrix, ())
    exact = DomainMatrix.from_Matrix(matrix, extension=True)
    exact = exact.convert_to(_algebraic_field(exact.domain))
    try:
        return exact.inv().to_Matrix()
    except DMNonInvertibleMatrixError as error:
        raise NonInvertibleMatrixError("Matrix is not invertible") from error


def normalize_polynomials(expressions, variables):
    """Return polynomials with canonical coefficients in a shared exact field.

    ``variables`` must list every genuine symbolic variable.  Coefficient
    identities, including algebraic zeros, are reduced exactly; powers of the
    listed variables are retained.  Empty lists and constant polynomials are
    supported.  Non-polynomial expressions are rejected by SymPy's polynomial
    conversion rather than silently changing their meaning.
    """
    variables = tuple(variables)
    if any(not isinstance(variable, sp.Symbol) for variable in variables):
        raise TypeError("Polynomial variables must be symbols")
    if len(set(variables)) != len(variables):
        raise ValueError("Polynomial variables must be distinct")
    values = _exact_expressions(expressions, variables)
    if not values:
        return []
    # SymPy requires at least one generator even for constant polynomials.
    generators = variables or (sp.Dummy("coefficient_generator"),)
    polynomials, options = sp.parallel_poly_from_expr(
        values, *generators, extension=True)
    field = _algebraic_field(options["domain"])
    return [polynomial.set_domain(field).as_expr() for polynomial in polynomials]
