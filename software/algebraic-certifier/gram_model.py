"""Exact affine Gram models from a closed body's combinatorics.

The origin is the vertex centroid, not the subsequently selected sphere center.
Convex metric uniqueness makes every full-body combinatorial automorphism an
orthogonal symmetry at this origin. All equations below are linear necessary
conditions on its Gram matrix; rank and positivity are deliberately separate.
Numerical positions supply only values of free parameters for branch selection.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
import time
import numpy as np
import sympy as sp

try:
    from .exact_linear import exact_inverse, normalize_polynomials
    from .patches import combinatorial_automorphisms, squared_polygon_distance
except ImportError:
    from exact_linear import exact_inverse, normalize_polynomials
    from patches import combinatorial_automorphisms, squared_polygon_distance


def _pair(i, j):
    return (min(i, j), max(i, j))


def _linear_reduce(equations, variables):
    """RREF over the exact algebraic coefficient field, with no heuristic zero."""
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    augmented = matrix.row_join(rhs)
    from sympy.polys.matrices import DomainMatrix
    dm = DomainMatrix.from_Matrix(augmented, extension=True).to_field()
    reduced, pivots = dm.rref()
    reduced = reduced.to_Matrix()
    if len(variables) in pivots:
        raise ValueError("Exact regular-face and symmetry constraints are inconsistent")
    free = [j for j in range(len(variables)) if j not in pivots]
    substitution = {}
    for row, column in enumerate(pivots):
        substitution[variables[column]] = sp.expand(
            reduced[row, -1] - sum(reduced[row, j] * variables[j] for j in free))
    return substitution, free, len(pivots), len(equations)


def build_gram_model(case, *, max_vertices=32, max_automorphisms=256,
                     max_flags=4096, include_face_affine=True):
    """Build G(t), satisfying centering, metric, symmetry and face flatness.

    ``case`` needs ``faces`` and may contain BI ``positions``. A returned
    ``witness`` is a numerical guide; no numerical quantity enters an equation.
    ``linear_basis`` contains the constant matrix followed by one coefficient
    matrix per returned variable. Labels need not be consecutive.
    """
    start = time.monotonic()
    faces = [list(map(int, f)) for f in case["faces"]]
    vertices = sorted(set(sum(faces, [])))
    n = len(vertices)
    if n > max_vertices:
        raise ValueError(f"Affine Gram vertex bound exceeded ({n}>{max_vertices})")
    index = {v: i for i, v in enumerate(vertices)}
    symmetry = combinatorial_automorphisms(
        faces, max_automorphisms=max_automorphisms, max_flags=max_flags)
    if not symmetry["permutations"]:
        raise ValueError("A valid closed regular-face sphere is required")
    permutations = [[index[v] for v in p] for p in symmetry["permutations"]]
    pairs = [_pair(i, j) for i in range(n) for j in range(i, n)]
    remaining = set(pairs)
    pair_orbits = []
    while remaining:
        seed = min(remaining)
        orbit, pending = {seed}, [seed]
        while pending:
            i, j = pending.pop()
            for permutation in permutations:
                image = _pair(permutation[i], permutation[j])
                if image not in orbit:
                    orbit.add(image)
                    pending.append(image)
        pair_orbits.append(sorted(orbit))
        remaining.difference_update(orbit)
    raw = sp.symbols(f"g0:{len(pair_orbits)}", real=True)
    at = {pair: symbol for orbit, symbol in zip(pair_orbits, raw) for pair in orbit}
    gram = sp.Matrix(n, n, lambda i, j: at[_pair(i, j)])
    equations = list(gram * sp.ones(n, 1))
    known, conflicts = {}, []
    for face in faces:
        m = len(face)
        for i, j in combinations(range(m), 2):
            pair = _pair(index[face[i]], index[face[j]])
            distance = squared_polygon_distance(m, j-i)
            if pair in known and sp.simplify(known[pair] - distance) != 0:
                conflicts.append((pair, known[pair], distance))
            known[pair] = distance
    if conflicts:
        raise ValueError("Conflicting exact distances within input faces")
    for (i, j), distance in known.items():
        equations.append(gram[i, i] + gram[j, j] - 2*gram[i, j] - distance)
    affine_relations = []
    if include_face_affine:
        for face in faces:
            m = len(face)
            if m <= 3:
                continue
            coefficient = sp.simplify(1 + 2*sp.cos(2*sp.pi/m))
            for k in range(m-3):
                vector = sp.zeros(n, 1)
                for offset, value in enumerate((-1, coefficient, -coefficient, 1)):
                    vector[index[face[k+offset]]] += value
                affine_relations.append(vector)
                equations.extend(list(gram*vector))
    equations = list(dict.fromkeys(sp.expand(e) for e in equations if e != 0))
    substitution, free_indices, rank, equation_count = _linear_reduce(equations, raw)
    free_raw = [raw[j] for j in free_indices]
    variables = sp.symbols(f"t0:{len(free_raw)}", real=True)
    rename = dict(zip(free_raw, variables))
    gram = gram.subs(substitution, simultaneous=True).subs(rename).applyfunc(sp.expand)
    zero = dict.fromkeys(variables, sp.S.Zero)
    constant = gram.subs(zero)
    basis = [constant] + [gram.diff(variable) for variable in variables]
    witness = {}
    if case.get("positions") is not None:
        points = np.asarray(case["positions"], dtype=float)
        if points.shape != (n, 3) or not np.isfinite(points).all():
            raise ValueError("Expected one finite numerical 3-vector per vertex")
        points = points - points.mean(axis=0)
        numerical_gram = points @ points.T
        for variable, free_index in zip(variables, free_indices):
            orbit = pair_orbits[free_index]
            witness[variable] = float(np.mean([numerical_gram[i, j] for i, j in orbit]))
    trace = [{"action": "affine_centroid_gram", "vertices": n,
              "automorphisms": len(permutations), "symmetry_complete": symmetry["complete"],
              "pair_orbits": len(pair_orbits), "linear_equations": equation_count,
              "linear_rank": rank, "free_parameters": len(variables),
              "regular_face_affine_relations": len(affine_relations),
              "exact_coefficients_only": True, "seconds": round(time.monotonic()-start, 6)}]
    return {"gram": gram, "variables": tuple(variables), "witness": witness,
            "linear_basis": basis, "known_constraints": known,
            "face_affine_relations": affine_relations,
            "vertices": vertices, "permutations": permutations,
            "symmetry": symmetry, "pair_orbits": pair_orbits,
            "parameter_orbits": [pair_orbits[j] for j in free_indices],
            "trace": trace}


def _compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def _independent_columns(matrix):
    from sympy.polys.matrices import DomainMatrix
    _, pivots = DomainMatrix.from_Matrix(matrix, extension=True).to_field().rref()
    return list(pivots)


def symmetry_blocks(model, *, positions=None, rank_tolerance=1e-7,
                    max_subgroups=512, max_generators=3):
    """Decompose G into rational character blocks of commuting involutions.

    Numerical rank is explicitly branch metadata, not a rigidity certificate.
    The exact independent verifier must accept any eventual solved branch.
    The character projectors and all matrix coefficients are exact.
    """
    n = model["gram"].rows
    identity = tuple(range(n))
    permutations = [tuple(p) for p in model["permutations"]]
    involutions = [p for p in permutations if p != identity and _compose(p,p)==identity]
    states = [((), frozenset([identity]))]
    seen = {frozenset([identity])}
    for generators, subgroup in states:
        if len(states) >= max_subgroups or len(generators) >= max_generators:
            continue
        for p in involutions:
            if p in subgroup or any(_compose(p,q) != _compose(q,p) for q in generators):
                continue
            expanded = subgroup | frozenset(_compose(p,q) for q in subgroup)
            if expanded in seen:
                continue
            if len(states) >= max_subgroups:
                break
            seen.add(expanded)
            states.append((generators+(p,), expanded))
    # More characters usually yield smaller rank equations. At equal order,
    # compare actual image ranks, using BI only as a candidate-choice guide.
    largest = max(len(subgroup) for _,subgroup in states)
    states = [(g,s) for g,s in states if len(s)==largest]
    points = None
    if positions is not None:
        points = np.asarray(positions,dtype=float).copy()
        points -= points.mean(axis=0)
    gram = model["gram"]
    candidates = []
    from itertools import product
    for generators, subgroup in states:
        matrices = [sp.Matrix(n,n,lambda i,j:int(p[j]==i)) for p in generators]
        blocks = []
        for signs in product((1,-1), repeat=len(generators)):
            projector = sp.eye(n)
            for sign, permutation in zip(signs,matrices):
                projector = projector*(sp.eye(n)+sign*permutation)/2
            columns = _independent_columns(projector)
            if not columns:
                continue
            basis = projector[:,columns]
            raw_block = (basis.T*gram*basis).applyfunc(sp.expand)
            variables = model["variables"]
            constant = raw_block.subs(dict.fromkeys(variables,0))
            coefficient_stack = sp.Matrix.vstack(constant,*(raw_block.diff(v) for v in variables))
            keep = _independent_columns(coefficient_stack)
            if not keep:
                continue
            basis = basis[:,keep]
            block = (basis.T*gram*basis).applyfunc(sp.expand)
            numerical_rank = None
            numerical_eigenvalues = None
            if points is not None:
                projected = np.asarray(basis,dtype=float).T @ points
                singular = np.linalg.svd(projected,compute_uv=False)
                numerical_rank = int(np.count_nonzero(singular > rank_tolerance))
                numerical_eigenvalues = (singular*singular).tolist()
            blocks.append({"character": list(signs), "basis": basis,
                           "gram": block, "dimension": block.rows,
                           "numerical_rank": numerical_rank,
                           "numerical_eigenvalues": numerical_eigenvalues,
                           "rank_is_branch_guide_only": True})
        # Prefer rank-one blocks and the smallest matrices; all such choices
        # describe the same exact full-body invariant Gram family.
        score = (sum((b["numerical_rank"] or 0)**2 for b in blocks),
                 sum(b["dimension"]**2 for b in blocks),tuple(generators))
        candidates.append((score,generators,blocks))
    _,generators,blocks = min(candidates,key=lambda item:item[0])
    model["blocks"] = blocks
    model["block_generators"] = [list(p) for p in generators]
    model["trace"].append({"action":"commuting_involution_character_blocks",
                           "subgroups_examined":len(states),
                           "subgroup_order":2**len(generators),
                           "block_dimensions":[b["dimension"] for b in blocks],
                           "rank_guides":[b["numerical_rank"] for b in blocks],
                           "coefficients_exact":True,
                           "numerical_ranks_are_not_proofs":True})
    return blocks


def anchor_parameterization(model, *, prefix="h"):
    """Choose exact affine coordinates that expose rank-one anchor products.

    Zero-rank block obligations are imposed linearly first. Independent anchor
    row entries are then made parameters, so Schur equations retain small
    products instead of products of large affine sums. This is a reversible
    exact change of variables on the selected zero-rank branch, not numerical
    coefficient recognition. The returned original_parameter_expressions map
    the input variables to the new ones.
    """
    blocks = model.get("blocks", [])
    original = tuple(model["variables"])
    zero_equations = [b["gram"][i,j] for b in blocks if b["numerical_rank"]==0
                      for i in range(b["dimension"]) for j in range(i,b["dimension"])]
    if zero_equations:
        zero_sub, free_indices, _, _ = _linear_reduce(zero_equations,original)
    else:
        zero_sub, free_indices = {}, list(range(len(original)))
    active = tuple(original[i] for i in free_indices)
    gram = model["gram"].xreplace(zero_sub).applyfunc(sp.expand)
    working_blocks=[]
    for block in blocks:
        working_blocks.append(dict(block,gram=block["gram"].xreplace(zero_sub).applyfunc(sp.expand)))
    candidates=[]
    anchors=[]
    for block_index,block in enumerate(working_blocks):
        if block["numerical_rank"] != 1:
            continue
        matrix=block["gram"]
        def anchor_score(i):
            diagonal=matrix[i,i]
            try:value=abs(float(diagonal.evalf(subs=model["witness"])))
            except (TypeError,ValueError):value=0.
            return (value<1e-10, bool(diagonal.free_symbols),sp.count_ops(diagonal),-value,i)
        anchor=min(range(matrix.rows),key=anchor_score)
        anchors.append({"block":block_index,"anchor":anchor})
        for j in range(matrix.cols):
            expression=matrix[anchor,j]
            if expression.free_symbols:
                candidates.append((int(bool(matrix[anchor,anchor].free_symbols)),expression,
                                   {"block":block_index,"row":anchor,"column":j,"kind":"anchor_entry"}))
    candidates.sort(key=lambda item:(item[0],sp.count_ops(item[1]),str(item[1])))
    for block_index,block in enumerate(working_blocks):
        matrix=block["gram"]
        for i in range(matrix.rows):
            for j in range(i,matrix.cols):
                if matrix[i,j].free_symbols:
                    candidates.append((2,matrix[i,j],{"block":block_index,"row":i,"column":j,"kind":"remaining_entry"}))
    candidates += [(3,v,{"kind":"previous_parameter","parameter":v}) for v in active]
    rows,forms,sources=[],[],[]
    for _,expression,source in candidates:
        row=[expression.diff(v) for v in active]
        if not any(row):
            continue
        trial=sp.Matrix(rows+[row])
        if len(_independent_columns(trial.T)) <= len(rows):
            continue
        rows.append(row);forms.append(expression);sources.append(source)
        if len(rows)==len(active):break
    variables=sp.symbols(f"{prefix}0:{len(active)}",real=True)
    if active:
        matrix=sp.Matrix(rows)
        offsets=sp.Matrix([e.subs(dict.fromkeys(active,0)) for e in forms])
        replacements=exact_inverse(matrix)*(sp.Matrix(variables)-offsets)
        parameter_sub=dict(zip(active,normalize_polynomials(replacements,variables)))
    else:
        parameter_sub={}
    substitution={v:sp.expand(zero_sub.get(v,v).xreplace(parameter_sub)) for v in original}
    result=dict(model)
    result["variables"]=tuple(variables)
    result["gram"]=gram.xreplace(parameter_sub).applyfunc(sp.expand)
    result["blocks"]=[dict(b,gram=b["gram"].xreplace(parameter_sub).applyfunc(sp.expand)) for b in working_blocks]
    result["witness"]={v:float(e.evalf(subs=model["witness"])) for v,e in zip(variables,forms)} if model["witness"] else {}
    result["linear_basis"]=[result["gram"].subs(dict.fromkeys(variables,0))]+[result["gram"].diff(v) for v in variables]
    result["original_parameter_expressions"]=substitution
    result["anchor_parameter_forms"]=forms
    result["trace"]=list(model["trace"])+[{"action":"exact_anchor_entry_parameterization",
        "rankzero_linear_equations":len(zero_equations),"old_parameters":len(original),
        "new_parameters":len(variables),"anchor_indices":anchors,"parameter_sources":sources,
        "parameter_expressions":substitution,"coefficients_exact":True}]
    return result


def rank_one_coordinate_system(model, *, max_coordinates=16):
    """Parameterize each rank-one character block as an outer product.

    Temporary coordinate signs are independent choices of symmetry axes. We
    ultimately retain only their quadratic products (the original Gram
    parameters), so those signs need not enlarge the recovered number field.
    """
    original=tuple(model["variables"])
    zero=dict.fromkeys(original,0)
    rows,rhs,variables=[],[],[]
    witness,axes={},[]
    for block_id,block in enumerate(model.get("blocks",[])):
        matrix=block["gram"]
        rank=block["numerical_rank"]
        dimension=matrix.rows
        if rank==0:
            target=sp.zeros(dimension)
            axes.append({"block":block_id,"rank":0})
        elif rank==1:
            constants=[i for i in range(dimension)
                       if not matrix[i,i].free_symbols and matrix[i,i].is_positive]
            coordinates=list(sp.symbols(f"axis{block_id}_0:{dimension}",real=True))
            if constants:
                anchor=constants[0]
                coordinates[anchor]=sp.sqrt(matrix[anchor,anchor])
            else:
                anchor=max(range(dimension),
                           key=lambda i:float(matrix[i,i].evalf(subs=model["witness"])))
            guide_norm=float(sp.sqrt(matrix[anchor,anchor]).evalf(subs=model["witness"]))
            if not np.isfinite(guide_norm) or guide_norm<=1e-10:
                raise ValueError("Rank-one block has no nonzero numerical anchor guide")
            for i,coordinate in enumerate(coordinates):
                if isinstance(coordinate,sp.Symbol):
                    variables.append(coordinate)
                    witness[coordinate]=float(matrix[anchor,i].evalf(subs=model["witness"]))/guide_norm
            vector=sp.Matrix(coordinates)
            target=vector*vector.T
            axes.append({"block":block_id,"rank":1,"anchor":anchor,
                         "coordinates":coordinates,"fixed_anchor":bool(constants),
                         "axis_sign_is_gauge":True})
        else:
            raise ValueError("Outer-product route requires only rank-zero/one character blocks")
        for i in range(dimension):
            for j in range(i,dimension):
                entry=matrix[i,j]
                rows.append([entry.diff(v) for v in original])
                rhs.append(sp.expand(target[i,j]-entry.subs(zero)))
    if len(variables)>max_coordinates:
        raise ValueError("Rank-one coordinate count exceeds its construction bound")
    linear=sp.Matrix(rows)
    target=sp.Matrix(rhs)
    pivot_rows=_independent_columns(linear.T)
    if len(pivot_rows)!=len(original):
        raise ValueError("Symmetry blocks do not determine every affine Gram parameter")
    if original:
        parameter_values=exact_inverse(linear[pivot_rows,:])*target[pivot_rows,:]
        parameter_values=sp.Matrix(normalize_polynomials(parameter_values,variables))
    else:
        parameter_values=sp.zeros(0,1)
    equations=[sp.expand(e) for e in target-linear*parameter_values]
    equations=list(dict.fromkeys(e for e in equations if e!=0))
    return {"equations":equations,"variables":tuple(variables),"witness":witness,
            "gram_parameter_expressions":dict(zip(original,map(sp.expand,parameter_values))),
            "trace":[{"action":"rank_one_outer_products","axes":axes,
                       "coordinate_variables":len(variables),"equations":len(equations),
                       "gram_parameters_eliminated_linearly":len(original),
                       "numerical_values_used_only_as_guides":True}]}


def polynomial_linear_span(equations,variables,witness,trace,*,max_rounds=10,max_monomials=512):
    """Expose affine consequences without increasing polynomial degrees.

    Row reduction regards the occurring monomials as formal independent
    columns. It only takes constant linear combinations of given equations;
    unlike an unrestricted rational substitution it cannot raise their degree.
    """
    try:
        from .reduce import _select_factors,_unique
    except ImportError:
        from reduce import _select_factors,_unique
    from sympy.polys.matrices import DomainMatrix
    variables=list(variables)
    equations=list(equations)
    substitutions={}
    for turn in range(max_rounds):
        equations=_unique([_select_factors(e,witness,trace) for e in equations])
        if not variables or not equations:
            return equations,variables,substitutions
        polynomials=[sp.Poly(e,*variables) for e in equations]
        monomials=sorted(set(m for p in polynomials for m,c in p.terms()),
                         key=lambda m:(sum(m),m),reverse=True)
        if len(monomials)>max_monomials:
            raise ValueError("Polynomial linear-span monomial bound exceeded")
        coefficients=sp.Matrix([[p.coeff_monomial(m) for m in monomials] for p in polynomials])
        reduced,pivots=DomainMatrix.from_Matrix(coefficients,extension=True).to_field().rref()
        reduced=reduced.to_Matrix()
        expressions=[sp.prod(v**power for v,power in zip(variables,m)) for m in monomials]
        equations=[sp.expand(sum(c*m for c,m in zip(reduced.row(i),expressions)))
                   for i in range(len(pivots))]
        equations=_unique([_select_factors(e,witness,trace) for e in equations])
        affine=[e for e in equations if sp.Poly(e,*variables).total_degree()==1]
        trace.append({"action":"polynomial_linear_span","round":turn,
                      "equations":len(equations),"variables":len(variables),
                      "affine_consequences":len(affine),"degree_increased":False})
        if not affine:
            return equations,variables,substitutions
        current,free,_,_=_linear_reduce(affine,variables)
        substitutions={**{v:sp.expand(e.xreplace(current)) for v,e in substitutions.items()},**current}
        equations=[sp.expand(e.xreplace(current)) for e in equations]
        variables=[variables[j] for j in free]
    return equations,variables,substitutions


def solve_rank_one_blocks(model,*,max_seconds=30,max_coordinates=16):
    """Return complete witnessed assignments of the original Gram parameters.

    This solves a small coordinate ideal, but recovers only invariant Gram
    products. Reducing those products modulo the selected irreducible factor
    removes temporary axis-sign fibers without numerical algebraic guessing.
    Every result still requires exact PSD/rank, center and body verification.
    """
    try:
        from .reduce import _time_limit,ReductionTimeout,_real_roots_near
    except ImportError:
        from reduce import _time_limit,ReductionTimeout,_real_roots_near
    started=time.monotonic()
    trace=[]
    result={"status":"UNRESOLVED","solutions":[],"derivation":trace}
    try:
        with _time_limit(max_seconds):
            system=rank_one_coordinate_system(model,max_coordinates=max_coordinates)
            trace.extend(system["trace"])
            equations,variables,substitutions=polynomial_linear_span(
                system["equations"],system["variables"],system["witness"],trace)
            if not variables:
                raise ValueError("No primitive coordinate remains for this bounded route")
            before=time.monotonic()
            basis=sp.groebner(equations,*variables,order="grevlex")
            trace.append({"action":"coordinate_ideal_groebner","variables":len(variables),
                          "equations":len(equations),"basis_size":len(basis.polys),
                          "input_polynomials":equations,"variable_symbols":variables,
                          "zero_dimensional":basis.is_zero_dimensional,
                          "seconds":time.monotonic()-before})
            if not basis.is_zero_dimensional:
                raise ValueError("Coordinate ideal is not zero-dimensional")
            lex=basis.fglm("lex")
            trace.append({"action":"coordinate_lex_basis","basis":list(lex.exprs)})
            univariates=[p.as_expr() for p in lex.polys if len(p.as_expr().free_symbols)==1]
            if not univariates:
                raise ValueError("Coordinate basis has no univariate polynomial")
            polynomial=min(univariates,key=lambda e:sp.Poly(e).degree())
            primitive=next(iter(polynomial.free_symbols))
            roots=_real_roots_near(polynomial,primitive,system["witness"][primitive])
            for root in roots[:4]:
                factor=(root.poly.as_expr().subs(root.poly.gen,primitive)
                        if isinstance(root,sp.CRootOf) else primitive-root)
                component=sp.groebner([p.as_expr() for p in lex.polys]+[factor],*variables,order="lex")
                expressions={}
                for variable,expression in system["gram_parameter_expressions"].items():
                    expression=sp.expand(expression.xreplace(substitutions))
                    expressions[variable]=component.reduce(expression)[1]
                if any(e.free_symbols-{primitive} for e in expressions.values()):
                    trace.append({"action":"metric_fiber_recovery","success":False,
                                  "reason":"a Gram parameter still varies in the chosen coordinate fiber"})
                    continue
                assignment={v:e.subs(primitive,root) for v,e in expressions.items()}
                if any(abs(float(e.evalf())-model["witness"][v])>1e-5*max(1,abs(model["witness"][v]))
                       for v,e in assignment.items()):
                    continue
                trace.append({"action":"metric_fiber_recovery","success":True,
                              "primitive_variable":primitive,"primitive_polynomial":factor,
                              "primitive_root":root,"metric_parameter_expressions":expressions,
                              "coordinate_axis_signs_removed":True})
                result["solutions"].append(assignment)
            if result["solutions"]:
                result["status"]="ALGEBRAIC_CANDIDATE"
            else:
                result["reason"]="No complete witnessed metric point recovered from the bounded coordinate ideal"
    except ReductionTimeout as exc:
        result.update(reason=str(exc),timed_out=True)
    except (ValueError,TypeError,NotImplementedError,ZeroDivisionError,
            sp.PolynomialError,sp.polys.polyerrors.CoercionFailed) as exc:
        result["reason"]=f"{type(exc).__name__}: {exc}"
    result["seconds"]=time.monotonic()-started
    return result
