"""Small exact Gram coordinates chosen from projected distance vectors."""
from itertools import combinations
import time
import sympy as sp

from gram_model import _linear_reduce, _independent_columns
from exact_linear import exact_inverse, normalize_polynomials
from patches import enrich_system
from reduce import _unique, _select_factors, _combine_substitutions


def _normalized_matrix(matrix, variables):
    return sp.Matrix(matrix.rows, matrix.cols,
                     normalize_polynomials(list(matrix), variables))


def add_reflection_patch_relations(model, case):
    """Add affine parallelogram identities from independently forced squares."""
    groups = case['groups']
    symbols = sp.symbols(f's0:{len(groups)}', positive=True)
    patch = enrich_system({'faces': case['faces'], 'groups': groups,
                           'radii_symbols': symbols, 'radius_variables': symbols,
                           'class_of': {v:i for i,g in enumerate(groups) for v in g}})
    equations = []
    index = {v:i for i,v in enumerate(model['vertices'])}
    for square in patch['reflection_squares']:
        vector = sp.zeros(model['gram'].rows, 1)
        for i,v in enumerate(square['vertices']):
            vector[index[v]] += (-1)**i
        equations.extend(list(model['gram']*vector))
    if not equations:
        return model
    sub, free, rank, _ = _linear_reduce(equations, model['variables'])
    result = dict(model)
    result['gram'] = model['gram'].xreplace(sub).applyfunc(sp.expand)
    result['variables'] = tuple(model['variables'][i] for i in free)
    result['linear_basis'] = [result['gram'].subs(dict.fromkeys(result['variables'],0))]
    result['linear_basis'] += [result['gram'].diff(v) for v in result['variables']]
    result['trace'] = list(model['trace'])+[{
        'operation':'forced_reflection_square_affine_relations',
        'squares':patch['reflection_squares'], 'independent_linear_constraints':rank,
        'parameter_substitutions':sub}]
    return result


def projected_distance_parameterization(model, *, prefix='k', max_pairs=256, max_anchor_trials=1024):
    """Use constant local metrics as Schur anchors whenever available.

    Candidate anchor vectors are all pair differences projected onto a symmetry
    character.  Constancy and independence of their metric are checked exactly;
    the BI guide cannot create a constant length or angle.  The remaining free
    affine forms are reparameterized by their dot products with these anchors.
    """
    n = model['gram'].rows
    gram, variables = model['gram'], tuple(model['variables'])
    zero = dict.fromkeys(variables, 0)
    permutations = [sp.Matrix(n,n,lambda i,j:int(p[j]==i))
                    for p in model['block_generators']]
    new_blocks=[]; sources=[]; forms=[]
    for block in model['blocks']:
        rank=block['numerical_rank']
        if rank not in (1,2,3):
            new_blocks.append(block)
            continue
        projector=sp.eye(n)
        for sign,p in zip(block['character'],permutations):
            projector=projector*(sp.eye(n)+sign*p)/2
        candidates=[]
        for pair_number,(i,j) in enumerate(combinations(range(n),2)):
            if pair_number>=max_pairs:break
            q=projector[:,i]-projector[:,j]
            if not any(q):continue
            norm=sp.expand((q.T*gram*q)[0])
            if norm.free_symbols or norm==0:continue
            if any(q==old[1] or q==-old[1] for old in candidates):continue
            candidates.append(((i,j),q,norm))
        chosen=None
        for trial_number,selection in enumerate(combinations(candidates,rank)):
            if trial_number>=max_anchor_trials:break
            B=sp.Matrix.hstack(*(item[1] for item in selection))
            metric=(B.T*gram*B).applyfunc(sp.expand)
            if any(value.free_symbols for value in metric):continue
            det=sp.simplify(metric.det())
            if det==0:continue
            score=(sum(sp.count_ops(e) for e in metric),
                   sum(sp.count_ops(e) for e in exact_inverse(metric)),
                   tuple(item[0] for item in selection))
            if chosen is None or score<chosen[0]:chosen=(score,selection,B,metric)
        if chosen is None:
            new_blocks.append(block)
            continue
        _,selection,B,metric=chosen
        # Extend these physical directions modulo common exact Gram nullspace.
        # This preserves the block family while putting the useful anchors first.
        directions=[B[:,i] for i in range(B.cols)]
        images=[]
        for direction in directions:
            column=gram*direction
            images.append(sp.Matrix.vstack(column.subs(zero),*(column.diff(v) for v in variables)))
        counts={i:sum(i in item[0] for item in selection) for i in range(n)}
        origin=min(range(n),key=lambda i:(-counts[i],i))
        extensions=[projector[:,i]-projector[:,origin] for i in range(n) if i!=origin]
        extensions += [block['basis'][:,i] for i in range(block['basis'].cols)]
        for direction in extensions:
            column=gram*direction
            image=sp.Matrix.vstack(column.subs(zero),*(column.diff(v) for v in variables))
            trial=sp.Matrix.hstack(*images,image)
            if len(_independent_columns(trial))>len(images):
                images.append(image);directions.append(direction)
        basis=sp.Matrix.hstack(*directions)
        matrix=(basis.T*gram*basis).applyfunc(sp.expand)
        updated=dict(block,basis=basis,gram=matrix,dimension=matrix.rows,
                     constant_anchor_indices=tuple(range(rank)))
        new_blocks.append(updated)
        sources.append({'character':block['character'],
                        'projected_pair_differences':[item[0] for item in selection],
                        'constant_metric':metric.tolist()})
        for i in range(rank):
            for j in range(matrix.cols):
                if matrix[i,j].free_symbols:forms.append(matrix[i,j])
    # Prefer sparse anchor products, then complete the exact affine coordinate basis.
    anchor_forms=set(forms)
    forms=sorted(anchor_forms,key=lambda e:(sp.count_ops(e),str(e)))
    for block in new_blocks:
        forms.extend(e for e in block['gram'] if e.free_symbols)
    forms.extend(variables)
    selected=[];rows=[]
    for expression in forms:
        row=[expression.diff(v) for v in variables]
        if not any(row):continue
        if len(_independent_columns(sp.Matrix(rows+[row]).T))>len(rows):
            selected.append(expression);rows.append(row)
        if len(rows)==len(variables):break
    new_variables=sp.symbols(f'{prefix}0:{len(variables)}',real=True)
    offset=sp.Matrix([e.subs(zero) for e in selected])
    vals=exact_inverse(sp.Matrix(rows))*(sp.Matrix(new_variables)-offset)
    sub=dict(zip(variables,normalize_polynomials(vals,new_variables)))
    result=dict(model)
    result['variables']=tuple(new_variables)
    result['gram']=_normalized_matrix(gram.xreplace(sub),new_variables)
    result['blocks']=[dict(b,gram=_normalized_matrix(b['gram'].xreplace(sub),new_variables)) for b in new_blocks]
    result['witness']={v:float(e.evalf(subs=model['witness'])) for v,e in zip(new_variables,selected)}
    result['original_parameter_expressions']=sub
    result['protected_projection_variables']=tuple(v for v,e in zip(new_variables,selected) if e in anchor_forms)
    result['trace']=list(model['trace'])+[{'operation':'constant_projected_distance_metric',
                                           'anchors':sources,'parameter_forms':selected,
                                           'parameter_substitutions':sub}]
    return result


def eliminate_origin_conic_pairs(equations, variables, witness, trace):
    """Solve small line/conic intersections that share the origin exactly.

    The line is homogeneous in a pair of coordinates. Substitution into a
    conic through the origin factors off the known zero solution and leaves a
    linear equation for the other intersection. The numerical guide chooses a
    candidate branch; all original equations and nonzero obligations are traced.
    """
    equations=list(equations);variables=list(variables);substitutions={}
    for _ in range(len(variables)//2):
        found=None
        for x,y in combinations(variables,2):
            pair={x,y}
            lines=[];conics=[]
            for equation in equations:
                if not pair<=equation.free_symbols:continue
                poly=sp.Poly(equation,x,y)
                if poly.coeff_monomial(1)!=0:continue
                if poly.total_degree()==1:lines.append(equation)
                if poly.total_degree()==2:conics.append(equation)
            for line in lines:
                for conic in conics:
                    if len((line.free_symbols|conic.free_symbols)-pair)>2:continue
                    coefficient=sp.diff(line,x)
                    try:
                        if abs(float(coefficient.evalf(subs=witness)))<1e-9:continue
                    except (ValueError,TypeError):continue
                    replacement=sp.cancel(-line.subs(x,0)/coefficient)
                    reduced=sp.factor(conic.subs(x,replacement))
                    factors=[_select_factors(reduced,witness,[])]
                    for factor in factors:
                        polynomial=sp.Poly(factor,y)
                        if polynomial.degree()!=1:continue
                        a,b=polynomial.all_coeffs()
                        yvalue=sp.cancel(-b/a)
                        xvalue=sp.cancel(replacement.subs(y,yvalue))
                        if xvalue.has(x,y) or yvalue.has(x,y):continue
                        score=(sp.count_ops(xvalue)+sp.count_ops(yvalue),str(x),str(y))
                        if found is None or score<found[0]:
                            found=(score,x,y,xvalue,yvalue,line,conic,coefficient,a)
        if found is None:break
        _,x,y,xvalue,yvalue,line,conic,coefficient,a=found
        step={x:xvalue,y:yvalue}
        trace.append({'operation':'line_conic_known_origin_branch',
                      'variables':[x,y],'line':line,'conic':conic,
                      'expressions':step,'nonzero_conditions':[coefficient,a],
                      'numerical_guide_selects_candidate_only':True})
        equations=_unique([sp.cancel(e.xreplace(step)).as_numer_denom()[0] for e in equations])
        variables=[v for v in variables if v not in step]
        substitutions=_combine_substitutions(substitutions,step)
    return equations,variables,substitutions



def eliminate_affine_equations(equations, variables, trace):
    """Finish exact degree-one constraints in their existing coefficient field.

    Projection variables are protected while nonlinear geometry is simplified,
    but a resulting affine system requires no root isolation. In particular,
    regular pentagon constants should remain in Q(sqrt(5)) rather than being
    sent through a rational polynomial lifting step.
    """
    variables=list(variables)
    if not variables:
        return list(equations),variables,{}
    affine=[e for e in equations if sp.Poly(e,*variables).total_degree()==1]
    if not affine:
        return list(equations),variables,{}
    substitutions,free,rank,_=_linear_reduce(affine,variables)
    reduced=_unique([sp.expand(e.xreplace(substitutions)) for e in equations])
    remaining=[variables[i] for i in free]
    trace.append({'operation':'exact_affine_projection_completion',
                  'equations':affine,'linear_rank':rank,
                  'expressions':substitutions,'remaining_variables':len(remaining),
                  'coefficient_field_preserved':True})
    return reduced,remaining,substitutions


def solve_projected_gram_model(case, *, max_seconds=45, max_vertices=32):
    """Bounded candidate solver using only generic local distance operations.

    This is useful when an invariant Gram block has rank greater than one.
    It discovers constant projected metrics and keeps their coordinates until
    local line/conic geometry has simplified the system. Exact verification
    remains the caller's responsibility.
    """
    from gram_model import build_gram_model, symmetry_blocks
    from distance_kernel import gram_rank_relations
    from reduce import _time_limit, ReductionTimeout, _eliminate_linear, _groebner_solutions
    from gram_solve import recover_radii
    started=time.monotonic();trace=[]
    result={'status':'UNRESOLVED','candidates':[],'derivation':trace}
    try:
        with _time_limit(max_seconds):
            model=add_reflection_patch_relations(build_gram_model(case,max_vertices=max_vertices),case)
            symmetry_blocks(model,positions=case['positions'])
            if sum(b['numerical_rank'] or 0 for b in model['blocks'])!=3:
                raise ValueError('No rank-three symmetry representation selected')
            model=projected_distance_parameterization(model)
            trace.extend(model['trace'])
            equations=[]
            for block in model['blocks']:
                relations=gram_rank_relations(block['gram'],block['numerical_rank'],witness=model['witness'])
                equations.extend(relations['constraints'])
                trace.append({'operation':'projected_local_gram_completion',
                              'rank':block['numerical_rank'],
                              'equations':relations['constraints'],
                              'nonzero_conditions':relations['nonzero_conditions'],
                              'positive_conditions':relations['positive_conditions']})
            variables=list(model['variables'])
            equations,variables,substitutions=_eliminate_linear(
                equations,variables,model['witness'],trace,
                protected=model['protected_projection_variables'])
            equations,variables,step=eliminate_origin_conic_pairs(
                equations,variables,model['witness'],trace)
            substitutions=_combine_substitutions(substitutions,step)
            equations,variables,step=eliminate_affine_equations(equations,variables,trace)
            substitutions=_combine_substitutions(substitutions,step)
            trace.append({'operation':'projected_reduced_system','variables':variables,
                          'equations':equations,'coordinate_parameters_preserved':True})
            solutions=_groebner_solutions(equations,variables,model['witness'],trace,max_variables=8)
            for solution in solutions:
                assignment=_combine_substitutions(substitutions,solution)
                if any(v not in assignment for v in model['variables']):continue
                result['candidates'].append(recover_radii(model,case,assignment,trace))
            if result['candidates']:
                result.update(status='ALGEBRAIC_CANDIDATE',
                              reason='Exact projected distance metrics and local line/conic completion')
            else:
                result['reason']='No isolated point found after generic projected distance completion'
    except ReductionTimeout as exc:
        result.update(reason=str(exc),timed_out=True)
    except (ValueError,TypeError,NotImplementedError,ZeroDivisionError,sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,sp.polys.polyerrors.NotInvertible) as exc:
        result['reason']=f'{type(exc).__name__}: {exc}'
    result['seconds']=time.monotonic()-started
    return result
