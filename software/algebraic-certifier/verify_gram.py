"""Independent exact certification of coordinates and arbitrary-center spheres.

Coordinates may use any exact basis with a positive definite metric.  The
certificate supplies data, never trusted claims: topology, all prescribed face
distances, supporting planes, and distances to the common sphere center are
recomputed here.  The center is allowed inside, on, or outside the polyhedron.
"""
from itertools import combinations
from collections import defaultdict

import sympy as sp
from sympy.matrices.exceptions import NonInvertibleMatrixError

try:
    from .verify import ExactField, ProofFailure, _topology, _face_coordinates
except ImportError:
    from verify import ExactField, ProofFailure, _topology, _face_coordinates


def _certify_facet_patch(face_ids, faces, points, coplanar_vertices, field):
    """Prove that convex face polygons tile one full supporting facet.

    Consistently oriented polygon boundaries must cancel to a single strictly
    convex boundary through exactly the vertices on the support plane.  The
    winding number of that boundary is the sum of the positive polygon winding
    numbers, proving coverage once and excluding gaps or overlapping interiors.
    _topology requires positive curvature, so valid input vertices cannot lie
    inside a maximal facet or in the relative interior of its boundary edges.
    """
    first = faces[face_ids[0]]
    p0, p1, p2 = (points[v] for v in first[:3])
    normal = (p1-p0).cross(p2-p0).applyfunc(field.normal)
    occurrences = defaultdict(list)
    for face_id in face_ids:
        face = list(faces[face_id])
        a, b, c = (points[v] for v in face[:3])
        orientation = field.sign(normal.dot((b-a).cross(c-a)))
        if orientation == 0:
            raise ProofFailure('A supporting facet contains a degenerate face')
        if orientation < 0:
            face.reverse()
        for i, a in enumerate(face):
            b = face[(i+1) % len(face)]
            occurrences[tuple(sorted((a,b)))].append((a,b))

    boundary_edges = []
    internal_edges = []
    for edge, directed in occurrences.items():
        if len(directed) == 1:
            boundary_edges.append(directed[0])
        elif len(directed) == 2 and directed[0] == directed[1][::-1]:
            internal_edges.append(list(edge))
        else:
            raise ProofFailure('Coplanar faces overlap or have inconsistent shared edges')
    outgoing, incoming = {}, {}
    for a,b in boundary_edges:
        if a in outgoing or b in incoming:
            raise ProofFailure('A coplanar facet boundary branches or repeats a vertex')
        outgoing[a] = b
        incoming[b] = a
    if not outgoing or set(outgoing) != set(incoming):
        raise ProofFailure('A coplanar facet has no closed boundary')
    start = min(outgoing)
    cycle, seen = [], set()
    vertex = start
    while vertex not in seen:
        seen.add(vertex)
        cycle.append(vertex)
        vertex = outgoing[vertex]
    if vertex != start or len(cycle) != len(boundary_edges):
        raise ProofFailure('A coplanar facet has multiple boundary components')
    if len(cycle) < 3 or set(cycle) != set(coplanar_vertices):
        raise ProofFailure('Coplanar faces do not cover the full supporting facet boundary')
    for i,a in enumerate(cycle):
        b = cycle[(i+1) % len(cycle)]
        for c in cycle:
            if c not in (a,b):
                turn = normal.dot((points[b]-points[a]).cross(points[c]-points[a]))
                if field.sign(turn) != 1:
                    raise ProofFailure('A coplanar facet boundary is not strictly convex')
    return {'faces':list(face_ids), 'vertices':list(coplanar_vertices),
            'boundary':cycle, 'internal_edges':sorted(internal_edges)}


def verify_gram_solution(faces, groups, squared_radii, geometric_certificate, *,
                         primitive_symbol=None, primitive_polynomial=None,
                         isolating_interval=None, primitive_interval=None):
    """Return CERTIFIED only after independent exact convex geometry checks.

    ``geometric_certificate`` contains ``metric`` (3 by 3), ``coordinates``
    (one three-vector per vertex), and ``center`` (one three-vector).  All values
    must be exact expressions in the specified real number field.  A numerical
    center, rank, face normal, or solver assertion cannot establish a check.

    Coplanar input faces are permitted only when their polygons provably tile
    an entire convex supporting facet.  The common sphere center may lie on a
    supporting plane or outside the convex body.
    """
    report={'status':'UNVERIFIED','verified':False,'certified':False,
            'method':'exact metric coordinates, supporting face planes, and arbitrary-center radii',
            'checks':{},'limitations':[]}
    try:
        if not isinstance(geometric_certificate,dict):
            raise ProofFailure('A geometric certificate object is required')
        if isolating_interval is None:
            isolating_interval=primitive_interval
        field=ExactField(primitive_symbol,primitive_polynomial,isolating_interval)

        def exact(value):
            value=sp.sympify(value)
            if value.has(sp.Float):
                raise ProofFailure('Geometric certificate entries must be exact, not floating point')
            return field.normal(value)

        metric=sp.Matrix(geometric_certificate['metric'])
        coordinates=sp.Matrix(geometric_certificate['coordinates'])
        center=sp.Matrix(geometric_certificate['center'])
        if metric.shape!=(3,3) or coordinates.cols!=3 or center.shape!=(3,1):
            raise ProofFailure('Expected a 3 by 3 metric, N by 3 coordinates, and a three-vector center')
        metric=metric.applyfunc(exact)
        coordinates=coordinates.applyfunc(exact)
        center=center.applyfunc(exact)
        if any(not field.zero(metric[i,j]-metric[j,i]) for i in range(3) for j in range(i)):
            raise ProofFailure('The coordinate metric is not symmetric')
        for k in range(1,4):
            field.positive(metric[:k,:k].det())
        report['checks']['positive_definite_coordinate_metric']=True
        n=coordinates.rows
        faces=[list(map(int,face)) for face in faces]
        groups=[list(map(int,group)) for group in groups]
        if len(groups)!=len(squared_radii):
            raise ProofFailure('The number of radii differs from the number of vertex groups')
        if any(not group for group in groups) or sorted(sum(groups,[]))!=list(range(n)):
            raise ProofFailure('Sphere groups must partition all coordinate vertices exactly once')
        edges=_topology(faces,n)
        report['checks']['combinatorial_sphere']=True
        report['checks']['positive_intrinsic_vertex_curvature']=n
        points=[coordinates.row(i).T for i in range(n)]

        def squared_norm(vector):
            return field.normal((vector.T*metric*vector)[0])

        distances={}
        for i,j in combinations(range(n),2):
            distance=squared_norm(points[i]-points[j])
            field.positive(distance)
            distances[i,j]=distance
        report['checks']['distinct_vertices']=n

        regular_pairs=0
        planes=[]
        center_signs=[]
        facet_groups=defaultdict(list)
        for face_id,face in enumerate(faces):
            reference=_face_coordinates(len(face))
            for i,j in combinations(range(len(face)),2):
                dx=reference[i][0]-reference[j][0]
                dy=reference[i][1]-reference[j][1]
                required=sp.simplify(dx*dx+dy*dy)
                pair=tuple(sorted((face[i],face[j])))
                if not field.zero(distances[pair]-required):
                    raise ProofFailure(f'Face {face_id} has an incorrect regular-face vertex distance')
                regular_pairs+=1
            p0,p1,p2=(points[v] for v in face[:3])
            u,v=p1-p0,p2-p0
            uv=field.normal((u.T*metric*v)[0])
            field.positive(squared_norm(u)*squared_norm(v)-uv*uv)
            covector=u.cross(v).applyfunc(field.normal)

            def side(point):
                return field.normal(covector.dot(point-p0))

            if any(not field.zero(side(points[vertex])) for vertex in face):
                raise ProofFailure(f'Face {face_id} is not planar')
            other=[vertex for vertex in range(n) if vertex not in face]
            if not other:
                raise ProofFailure('A face contains all vertices; the body is not three dimensional')
            signs=[field.sign(side(points[vertex])) for vertex in other]
            nonzero=[sign for sign in signs if sign!=0]
            if not nonzero:
                raise ProofFailure('All vertices lie on a face plane; the body is not three dimensional')
            inward=nonzero[0]
            if any(sign!=inward for sign in nonzero):
                raise ProofFailure(f'Face {face_id} does not support the convex hull')
            coplanar=sorted(face+[vertex for vertex,sign in zip(other,signs) if sign==0])
            facet_groups[tuple(coplanar)].append(face_id)
            center_sign=inward*field.sign(side(center))
            center_signs.append(center_sign)
            planes.append({'face':face_id,'inward_covector':[inward*x for x in covector],
                           'inward_offset':field.normal(inward*covector.dot(p0)),
                           'strictly_inner_vertices':len(nonzero),
                           'coplanar_vertices':coplanar,
                           'center_side':center_sign})
        facets=[_certify_facet_patch(ids,faces,points,vertices,field)
                for vertices,ids in facet_groups.items()]
        report['checks']['regular_face_pair_distances']=regular_pairs
        report['checks']['planar_regular_faces']=len(faces)
        report['checks']['supporting_face_planes']=len(faces)
        report['checks']['strictly_supporting_face_planes']=sum(
            len(plane['coplanar_vertices'])==len(faces[plane['face']]) for plane in planes)
        report['checks']['maximal_supporting_facets']=len(facets)
        report['checks']['coplanar_subdivided_facets']=sum(len(facet['faces'])>1 for facet in facets)
        report['checks']['exact_facet_tilings']=len(facets)
        report['checks']['embedded_closed_convex_boundary']=True

        radii=list(map(exact,squared_radii))
        for group,radius in zip(groups,radii):
            field.positive(radius)
            for vertex in group:
                if not field.zero(squared_norm(points[vertex]-center)-radius):
                    raise ProofFailure(f'Vertex {vertex} does not belong to its proposed sphere')
        report['checks']['positive_squared_radii']=True
        report['checks']['exact_sphere_memberships']=n
        location=('exterior' if -1 in center_signs else
                  'boundary' if 0 in center_signs else 'interior')
        report['center_location']=location
        report['center_face_signs']=center_signs
        report['supporting_faces']=planes
        report['supporting_facets']=facets
        report['status']='CERTIFIED'
        report['verified']=True
        report['certified']=True
        report['conclusion']=('The supplied exact coordinates form the prescribed closed convex '
                              'polyhedron with regular unit-edge faces, and its vertices lie '
                              f'on the proposed concentric spheres. The common center is {location}.')
        if field.interval is not None:
            report['refined_primitive_interval']=field.interval.strings()
    except (ProofFailure,ValueError,TypeError,KeyError,IndexError,NotImplementedError,
            ZeroDivisionError,sp.PolynomialError,NonInvertibleMatrixError) as exc:
        report['reason']=str(exc)
    return report
