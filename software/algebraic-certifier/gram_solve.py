"""Generic exact distance completion in symmetry-adapted Gram coordinates.

No body names, coordinate templates, or answer polynomials are used. Numerical
positions choose symmetry-rank branches and pivots; exact distances supply all
coefficients. Every returned radius candidate requires the independent verifier.
"""
from itertools import combinations
import time

import numpy as np
import sympy as sp

from gram_model import build_gram_model, symmetry_blocks, anchor_parameterization, solve_rank_one_blocks
from distance_kernel import gram_rank_relations
from reduce import (_time_limit, ReductionTimeout, _unique,
                    _primitive_candidate, solve_low_degree)
from verify import ExactField


def _field_rref(rows, columns, field):
    rows=[[field.normal(x) for x in row] for row in rows]
    pivots=[]
    for column in range(columns):
        pivot=next((i for i in range(len(pivots),len(rows)) if rows[i][column]!=0),None)
        if pivot is None:
            continue
        at=len(pivots)
        rows[at],rows[pivot]=rows[pivot],rows[at]
        divisor=rows[at][column]
        rows[at]=[field.normal(x/divisor) for x in rows[at]]
        for i in range(len(rows)):
            if i==at or rows[i][column]==0:
                continue
            multiplier=rows[i][column]
            rows[i]=[field.normal(x-multiplier*y) for x,y in zip(rows[i],rows[at])]
        pivots.append(column)
    if any(all(row[j]==0 for j in range(columns)) and row[-1]!=0 for row in rows):
        raise ValueError('Requested sphere partition is inconsistent with the exact completed metric')
    return rows,pivots


def recover_radii(model, case, assignment, trace):
    """Choose the center by affine equal-distance equations in the exact field."""
    variables=model['variables']
    primitive=_primitive_candidate([assignment[v] for v in variables],variables,{},trace)
    field=ExactField(primitive['primitive_symbol'],primitive['primitive_polynomial'],
                     primitive['isolating_interval'])
    replacements=primitive['radius_expressions']
    G=model['gram'].applyfunc(lambda e:field.normal(e.xreplace(replacements)))
    points=np.asarray(case['positions'],dtype=float)
    points=points-points.mean(axis=0)
    guide=points@points.T
    choices=sorted(combinations(range(len(points)),3),
                   key=lambda I:-abs(np.linalg.det(guide[np.ix_(I,I)])))
    anchors=None
    for I in choices:
        A=G.extract(I,I)
        if field.normal(A.det())!=0:
            anchors=I
            break
    if anchors is None:
        raise ValueError('Exact completed metric does not span three dimensions')
    # O=sum_j c_j p_anchor_j. Differences of equal squared radii are linear in c.
    rows=[]
    for group in case['groups']:
        base=group[0]
        for vertex in group[1:]:
            rows.append([2*(G[vertex,j]-G[base,j]) for j in anchors]
                        +[G[vertex,vertex]-G[base,base]])
    reduced,pivots=_field_rref(rows,3,field)
    free=[j for j in range(3) if j not in pivots]
    c0=sp.zeros(3,1)
    N=sp.zeros(3,len(free))
    for j,k in enumerate(free): N[k,j]=1
    for row,pivot in enumerate(pivots):
        c0[pivot]=reduced[row][-1]
        for j,k in enumerate(free): N[pivot,j]=-reduced[row][k]
    A=G.extract(anchors,anchors)
    coefficients=c0
    if free:
        # The nearest center to the vertex centroid is a deterministic exact
        # choice among the remaining affine family. A geometric certificate
        # also permits centers on or outside the realized body.
        normal=(N.T*A*N).applyfunc(field.normal)
        rhs=(-N.T*A*c0).applyfunc(field.normal)
        rr,pp=_field_rref([list(normal.row(i))+[rhs[i]] for i in range(len(free))],len(free),field)
        if len(pp)!=len(free):
            raise ValueError('Center minimum-norm system is singular')
        t=sp.Matrix([rr[i][-1] for i in range(len(free))])
        coefficients=(c0+N*t).applyfunc(field.normal)
    norm=field.normal((coefficients.T*A*coefficients)[0])
    inner=[field.normal(sum(G[i,j]*coefficients[k] for k,j in enumerate(anchors)))
           for i in range(len(points))]
    all_radii=[field.normal(G[i,i]-2*inner[i]+norm) for i in range(len(points))]
    radii=[]
    for group in case['groups']:
        radius=all_radii[group[0]]
        if any(field.normal(all_radii[v]-radius)!=0 for v in group):
            raise ValueError('Exact center failed a requested radius equality')
        if field.sign(radius)!=1:
            raise ValueError('A completed radius is not positive')
        radii.append(radius)
    # Check exact PSD/rank via one positive 3x3 principal block and its Schur
    # complement. This validates the numerical representation-rank selection.
    for k in range(1,4):
        if field.sign(A[:k,:k].det())!=1:
            raise ValueError('Completed Gram anchor is not positive definite')
    adj=A.adjugate().applyfunc(field.normal)
    determinant=field.normal(A.det())
    for i in range(G.rows):
        for j in range(i,G.cols):
            scalar=field.normal((G.extract([i],anchors)*adj*G.extract(anchors,[j]))[0])
            if field.normal(G[i,j]*determinant-scalar)!=0:
                raise ValueError('Completed Gram matrix failed exact rank-three verification')
    # Coordinates in the independent anchor-vector basis. Its positive metric
    # is A; no coordinate square roots or numerical recognition are needed.
    coordinates=[]
    for i in range(G.rows):
        vector=adj*G.extract(anchors,[i])
        coordinates.append([field.normal(x/determinant) for x in vector])
    geometry=dict(metric=A.tolist(), coordinates=coordinates,
                  center=list(coefficients))
    trace.append(dict(operation='exact_affine_sphere_center',anchors=anchors,
                      affine_center_rank=len(pivots),free_center_dimensions=len(free),
                      center_coefficients=list(coefficients),center_squared_norm=norm,
                      gram_rank_three_and_positive_semidefinite=True))
    symbols=sp.symbols(f's0:{len(radii)}',positive=True)
    return dict(status='ALGEBRAIC_CANDIDATE',primitive_symbol=primitive['primitive_symbol'],
                primitive_polynomial=primitive['primitive_polynomial'],
                isolating_interval=primitive['isolating_interval'],
                squared_radii=radii,radius_expressions=dict(zip(symbols,radii)),
                derivation=list(trace),auxiliary_assignments={},
                geometric_certificate=geometry)


def solve_gram_model(case, *, max_seconds=45, max_vertices=32):
    started=time.monotonic()
    trace=[]
    result=dict(status='UNRESOLVED',candidates=[],derivation=trace)
    try:
        with _time_limit(max_seconds):
            model=build_gram_model(case,max_vertices=max_vertices)
            symmetry_blocks(model,positions=case['positions'])
            if any((block['numerical_rank'] or 0)>1 for block in model['blocks']):
                from gram_reduction import solve_projected_gram_model
                return solve_projected_gram_model(case,
                    max_seconds=max(.1,max_seconds-(time.monotonic()-started)),
                    max_vertices=max_vertices)
            model=anchor_parameterization(model)
            blocks=model['blocks']
            trace.extend(model['trace'])
            ranks=[b['numerical_rank'] for b in blocks]
            if any(r is None for r in ranks) or sum(ranks)!=3:
                raise ValueError('No supported rank-three symmetry representation selected')
            witness=model['witness']
            equations=[]
            conditions=[]
            for index,block in enumerate(blocks):
                completion=gram_rank_relations(block['gram'],block['numerical_rank'],witness=witness)
                equations.extend(completion['constraints'])
                conditions.extend(completion['nonzero_conditions'])
                trace.append(dict(operation='local_gram_completion',block=index,
                    matrix=block['gram'].tolist(),rank_branch=block['numerical_rank'],
                    equations=completion['constraints'],conditions=completion['nonzero_conditions'],
                    anchors=completion['anchor_indices'],provenance=completion['provenance']))
            equations=_unique(equations)
            variables=list(model['variables'])
            if len(variables)>18:
                raise ValueError('Affine Gram parameter bound exceeded')
            remaining_budget=max(.1,max_seconds-(time.monotonic()-started))
            if len(variables)>7 and all(b['numerical_rank'] in (0,1) for b in blocks):
                reduction=solve_rank_one_blocks(model,max_seconds=remaining_budget)
            else:
                reduction=solve_low_degree(equations,variables,witness,
                                           max_seconds=remaining_budget,max_variables=8)
            trace.extend(reduction['derivation'])
            for assignment in reduction['solutions']:
                if any(v not in assignment for v in model['variables']):
                    continue
                candidate=recover_radii(model,case,assignment,trace)
                result['candidates'].append(candidate)
            if result['candidates']:
                result.update(status='ALGEBRAIC_CANDIDATE',reason='Generic exact Gram completion and affine center recovery')
            else:
                result['reason']='No isolated metric point extracted from the selected symmetry rank branch'
    except ReductionTimeout as exc:
        result['reason']=str(exc)
        result['timed_out']=True
    except (ValueError,TypeError,NotImplementedError,ZeroDivisionError,sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,sp.polys.polyerrors.NotInvertible) as exc:
        result['reason']=f'{type(exc).__name__}: {exc}'
    result['seconds']=time.monotonic()-started
    return result
