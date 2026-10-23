"""Numerical counterpart of the exact pigeonhole sphere search.

This searches equality-plane candidates using tolerances. It does NOT certify
the global optimum for thick spherical shells or for the underlying exact CRP.
Every returned grouping is rechecked against the requested radius diameter.
"""
import itertools
import time
import numpy as np


def _groups(points, center, tolerance):
    distances = np.linalg.norm(points-center, axis=1)
    order = np.argsort(distances)
    groups = []
    for index in order:
        i = int(index)
        if not groups or distances[i]-distances[groups[-1][0]] > tolerance:
            groups.append([i])
        else:
            groups[-1].append(i)
    radii = [float(np.mean(distances[g])) for g in groups]
    spreads = [float(np.ptp(distances[g])) for g in groups]
    return groups, radii, spreads


def solve_min_spheres_numeric(vertices, *, tolerance=1e-6, centers=(), max_nodes=20000, time_limit=2.):
    points = np.asarray(vertices, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points)==0 or not np.isfinite(points).all():
        raise ValueError('Expected finite 3D coordinates')
    if tolerance<=0 or not np.isfinite(tolerance):
        raise ValueError('Tolerance must be positive and finite')
    start = time.monotonic()
    shift = points.mean(axis=0)
    scale = float(np.max(np.linalg.norm(points-shift,axis=1))) or 1.
    P = (points-shift)/scale
    norms = np.einsum('ij,ij->i',P,P)
    tol = tolerance/scale
    # Rank tolerance is separate from the requested radius spread.
    rank_tol = max(1e-11,tol*.01)
    identity_tol = max(1e-10,tol*.05)
    nodes, skipped, completed_through = 0, 0, 0
    best = None
    decisions = []

    def score(c, refine=False):
        nonlocal best
        if not np.isfinite(c).all():
            return len(P)+1
        groups,radii,spreads = _groups(P,c,tol)
        if refine and any(len(g)>1 for g in groups):
            A,b = [],[]
            for g in groups:
                if len(g)>1:
                    p=P[g];s=norms[g]
                    A.extend(2*(p-p.mean(axis=0)))
                    b.extend(s-s.mean())
            A,b=np.asarray(A),np.asarray(b)
            step=np.linalg.lstsq(A,b-A@c,rcond=rank_tol)[0]
            trial=c+step
            other=_groups(P,trial,tol)
            if (len(other[0]),max(other[2])) <= (len(groups),max(spreads)):
                c=trial;groups,radii,spreads=other
        key=(len(groups),max(spreads),float(np.dot(c,c)))
        if best is None or key<best['key']:
            best=dict(key=key,center=c.copy(),groups=groups,radii=radii,spreads=spreads)
        return len(groups)

    score(np.zeros(3),True)
    for c in centers:
        score((np.asarray(c,dtype=float)-shift)/scale,True)
    # A direct linear fit handles the important cospherical case immediately.
    c=np.linalg.lstsq(2*(P-P.mean(axis=0)),norms-norms.mean(),rcond=rank_tol)[0]
    score(c,True)
    budget_hit=False

    def visit(origin,basis,k):
        nonlocal nodes,skipped
        if nodes>=max_nodes or time.monotonic()-start>time_limit:
            raise TimeoutError
        nodes+=1
        count=score(origin)
        if count<=k:
            score(origin,True)
            return True
        if basis.shape[1]==0:
            return False
        linear=-2*(P@basis)
        constant=norms-2*(P@origin)
        # Use a conservative numerical identity threshold. Final radii use tol.
        constant_tol=identity_tol*max(1.,float(np.linalg.norm(origin)))
        representatives=[]
        for i in range(len(P)):
            if any(np.linalg.norm(linear[i]-linear[j])<=identity_tol and abs(constant[i]-constant[j])<=constant_tol for j in representatives):
                continue
            representatives.append(i)
            if len(representatives)>k:
                break
        if len(representatives)<=k:
            # Approximate function identities alone do not override the actual
            # radius-spread check already performed above.
            skipped+=1
            return False
        for i,j in itertools.combinations(representatives,2):
            a=2*(P[j]-P[i])@basis
            rhs=norms[j]-norms[i]-2*np.dot(P[j]-P[i],origin)
            length=float(np.linalg.norm(a))
            if length<=rank_tol:
                skipped+=1
                continue
            direction=a/length
            child=origin+basis@direction*(rhs/length)
            null=np.linalg.svd(direction[None,:],full_matrices=True)[2][1:].T
            if visit(child,basis@null,k):
                return True
        return False

    for k in range(1,best['key'][0]):
        before=nodes
        try:
            found=visit(np.zeros(3),np.eye(3),k)
        except TimeoutError:
            budget_hit=True
            decisions.append(dict(k=k,outcome='budget',nodes=nodes-before))
            break
        decisions.append(dict(k=k,outcome='witness' if found else 'no_candidate',nodes=nodes-before))
        if found:
            break
        completed_through=k
        if best['key'][0]<=k+1:
            break

    score(best['center'],True)
    center=best['center']*scale+shift
    groups,radii,spreads=_groups(points,center,tolerance)
    if max(spreads)>tolerance*(1+1e-8):
        raise ArithmeticError('Returned groups fail the requested radius spread')
    return dict(k=len(groups),center=center.tolist(),radii=radii,
                squared_radii=[r*r for r in radii],groups=groups,radius_spreads=spreads,
                max_radius_spread=max(spreads),radius_tolerance=tolerance,
                minimum_certified=False,numerical_candidate_search_finished=not budget_hit,
                completed_candidate_tests_through=completed_through,
                approximate_or_parallel_branches_skipped=skipped,
                nodes=nodes,seconds=time.monotonic()-start,decisions=decisions,
                status='BUDGET_LIMIT' if budget_hit else 'NUMERICAL_CANDIDATE_RESULT',
                note='Smallest numerical sphere count found; not an exact minimality certificate or a complete thick-shell optimization.')
