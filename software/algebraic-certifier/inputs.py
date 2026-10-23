"""Read census/BI files. Numerical data are guides, never exact coefficients."""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def read_json(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_faces(faces):
    faces = [list(map(int, f)) for f in faces]
    if not faces or any(len(f) not in (3, 4, 5, 6, 8, 10) for f in faces):
        raise ValueError("Expected regular faces with 3, 4, 5, 6, 8, or 10 vertices")
    if any(len(set(f)) != len(f) for f in faces):
        raise ValueError("A face repeats a vertex")
    vertices = sorted(set(sum(faces, [])))
    if vertices != list(range(len(vertices))):
        raise ValueError("Vertex labels must be consecutive integers starting at zero")
    edges = Counter(tuple(sorted((f[j], f[(j+1) % len(f)])))
                    for f in faces for j in range(len(f)))
    if any(count != 2 for count in edges.values()):
        raise ValueError("Every edge of the closed input must belong to two faces")
    if len(vertices) - len(edges) + len(faces) != 2:
        raise ValueError("Input must have sphere Euler characteristic")
    for i in vertices:
        link = defaultdict(list)
        for f in faces:
            if i in f:
                at = f.index(i)
                j, k = f[at-1], f[(at+1) % len(f)]
                link[j].append(k)
                link[k].append(j)
        if not link or any(len(neighbors) != 2 for neighbors in link.values()):
            raise ValueError(f"Vertex {i} does not have a circular link")
        reached, pending = set(), [next(iter(link))]
        while pending:
            v = pending.pop()
            if v not in reached:
                reached.add(v)
                pending.extend(link[v])
        if len(reached) != len(link):
            raise ValueError(f"Vertex {i} has a disconnected link")
    return faces


def group_at_center(points, center, tolerance):
    distances = np.linalg.norm(points-center, axis=1)
    groups = []
    for i in np.argsort(distances):
        i = int(i)
        if not groups or distances[i]-distances[groups[-1][0]] > tolerance:
            groups.append([i])
        else:
            groups[-1].append(i)
    return dict(center=np.asarray(center).tolist(), groups=groups,
                squared_radii=[float(np.mean(distances[g]**2)) for g in groups],
                radius_spreads=[float(np.ptp(distances[g])) for g in groups])


def interior_clearance(points, faces, center):
    """Numerical preflight only. The exact verifier has final authority."""
    reference = points.mean(axis=0)
    clearances = []
    for face in faces:
        p, q, r = points[face[:3]]
        normal = np.cross(q-p, r-p)
        norm = np.linalg.norm(normal)
        if norm < 1e-12:
            return -float('inf')
        normal /= norm
        if np.dot(reference-p, normal) > 0:
            normal = -normal
        clearances.append(-float(np.dot(center-p, normal)))
    return min(clearances)


def load_case(case, *, bi=None, spheres=None, tolerance=1e-6):
    path = Path(case).expanduser().resolve()
    if path.is_dir():
        face_path, bi_path = path/'input.json', path/'numerical.json'
    else:
        face_path, bi_path = path, path.parent/'numerical.json'
    if bi:
        bi_path = Path(bi).expanduser().resolve()
    data = read_json(face_path)
    if isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict):
        data = data[0]
    faces = validate_faces(data['faces'] if isinstance(data, dict) else data)
    numeric = read_json(bi_path)
    if isinstance(numeric, list):
        numeric = next((r for r in numeric if r.get('success') and 'positions' in r),
                       next((r for r in numeric if 'positions' in r), {}))
    points = np.asarray(numeric.get('positions', numeric.get('vertices')), dtype=float)
    if points.shape != (1+max(sum(faces, [])), 3) or not np.isfinite(points).all():
        raise ValueError("BI file must contain finite positions for every input vertex")
    notes = []
    sources = {"faces": {"path": str(face_path), "sha256": digest(face_path)},
               "bi": {"path": str(bi_path), "sha256": digest(bi_path)}}
    if spheres:
        sphere_path = Path(spheres).expanduser().resolve()
        partition = read_json(sphere_path)
        sources['partition'] = {'path': str(sphere_path), 'sha256': digest(sphere_path)}
        if partition.get('input_sha256') and partition['input_sha256'] != digest(bi_path):
            raise ValueError("Sphere partition was computed from a different BI file (SHA256 mismatch)")
    else:
        from numeric_partition import solve_min_spheres_numeric
        seeds = [numeric['apex']] if 'apex' in numeric else []
        partition = solve_min_spheres_numeric(points, tolerance=tolerance,
                                               centers=seeds, time_limit=2.)
        if interior_clearance(points, faces, np.asarray(partition['center'])) <= 1e-8:
            notes.append("Sphere search selected a noninterior center; retaining its partition for exact metric verification")
    groups = [sorted(map(int, g)) for g in partition['groups']]
    if any(not g for g in groups) or sorted(sum(groups, [])) != list(range(len(points))):
        raise ValueError("Sphere groups must partition all vertices exactly once")
    center = np.asarray(partition['center'], dtype=float)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("Sphere center must be a finite 3-vector")
    distances = np.linalg.norm(points-center, axis=1)
    spread = max(float(np.ptp(distances[g])) for g in groups)
    if spread > tolerance:
        raise ValueError(f"Sphere partition does not match BI coordinates: radius spread {spread:g}")
    radius_values = [float(np.mean(distances[g]**2)) for g in groups]
    clearance = interior_clearance(points, faces, center)
    pair_error = 0.
    for f in faces:
        n = len(f)
        for j in range(n):
            for k in range(j):
                expected = np.sin(np.pi*(j-k)/n)/np.sin(np.pi/n)
                pair_error = max(pair_error, abs(np.linalg.norm(points[f[j]]-points[f[k]])-expected))
    differences = np.asarray([points[i]-points[g[0]] for g in groups for i in g[1:]])
    center_rank = int(np.linalg.matrix_rank(differences, tol=1e-7)) if len(differences) else 0
    return dict(faces=faces, groups=groups, squared_radii_numeric=radius_values,
                positions=points.tolist(), center=center.tolist(), sources=sources,
                name=partition.get('name', data.get('name', face_path.parent.name) if isinstance(data, dict) else face_path.parent.name),
                preflight={'numerical_only': True, 'center_clearance': clearance,
                           'radius_spread': spread, 'face_pair_error': float(pair_error),
                           'within_shell_center_rank': center_rank}, notes=notes)
