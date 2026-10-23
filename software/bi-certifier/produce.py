"""BI-seeded numerical Newton proposals. No coordinate recognition or symmetry routes."""

import json
import shutil
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import metric as m
import exact_intervals as iv
from check import precision, verify


def mid(x):
    return float(Q(x.lo + x.hi, 2 * iv.SCALE))


def propose_root(polygons, numerical, digits):
    cert = dict(
        method="interval-root",
        bits=max(160, 4 * digits + 80),
        polygons=polygons,
        flips=numerical["flips"],
        triangles=numerical["triangles"],
        base=numerical["triangles"][0],
    )
    precision(cert["bits"])
    V, faces, L, fixed, unknown, equations = m.setup(cert)
    p = numerical_frame(numerical)
    for v, row in fixed.items():
        p[v] = [mid(x) for x in row]
    index = {vk: i for i, vk in enumerate(unknown)}
    target = np.array([mid(L[e].square()) for e in equations])

    def fj(p):
        f = []
        J = np.zeros((len(equations), len(unknown)))
        for row, (v, w) in enumerate(equations):
            diff = p[v] - p[w]
            f.append(np.dot(diff, diff) - target[row])
            for vertex, sign in ((v, 1), (w, -1)):
                for k in range(3):
                    if (vertex, k) in index:
                        J[row, index[vertex, k]] = 2 * sign * diff[k]
        return np.array(f), J

    for _ in range(8):
        f, J = fj(p)
        step = np.linalg.solve(J, -f)
        for i, (v, k) in enumerate(unknown):
            p[v, k] += step[i]
        if max(abs(step)) < 1e-15:
            break
    _, J = fj(p)
    Y = np.linalg.inv(J)
    cert["preconditioner"] = [[str(Q(float(x))) for x in row] for row in Y]
    center = [str(Q(float(p[v, k]))) for v, k in unknown]
    if digits > 11:
        import mpmath as mp

        mp.mp.dps = digits + 35
        cv = lambda x: mp.mpf(x.lo + x.hi) / (2 * iv.SCALE)
        pm = [[mp.mpf(str(x)) for x in row] for row in p]
        for v, row in fixed.items():
            pm[v] = [cv(x) for x in row]
        targetm = [cv(L[e].square()) for e in equations]
        for _ in range(6):
            jm = mp.matrix(len(equations))
            fm = mp.matrix(len(equations), 1)
            for row, (v, w) in enumerate(equations):
                diff = [pm[v][k] - pm[w][k] for k in range(3)]
                fm[row] = sum(x * x for x in diff) - targetm[row]
                for vertex, sign in ((v, 1), (w, -1)):
                    for k in range(3):
                        if (vertex, k) in index:
                            jm[row, index[vertex, k]] = 2 * sign * diff[k]
            step = mp.lu_solve(jm, -fm)
            for i, (v, k) in enumerate(unknown):
                pm[v][k] += step[i]
            if max(abs(x) for x in step) < mp.mpf(10) ** (-digits - 15):
                break
        center = [mp.nstr(pm[v][k], digits + 25) for v, k in unknown]
    cert["center"] = center
    cert["radius"] = str(Q(1, 10 ** (digits + 1)))
    return cert, p


def numerical_frame(numerical):
    p = np.array(numerical["positions"], dtype=float)
    a, b, c = numerical["triangles"][0]
    origin = p[a].copy()
    e1 = p[b] - origin
    e1 /= np.linalg.norm(e1)
    e2 = p[c] - origin
    e2 -= e1 * np.dot(e1, e2)
    e2 /= np.linalg.norm(e2)
    p = (p - origin) @ np.array([e1, e2, np.cross(e1, e2)]).T
    if not np.all(np.isfinite(p)):
        raise ValueError("BI numerical gauge is degenerate")
    return p


def propose_metric(directory, polygons, numerical, digits):
    if not numerical.get("initialized") or "positions" not in numerical:
        raise ValueError("BI did not provide numerical coordinates")
    p = numerical_frame(numerical)
    report = None
    errors = []
    try:
        cert, p = propose_root(polygons, numerical, digits)
        (directory / "root-certificate.json").write_text(json.dumps(cert) + "\n")
        report = verify(cert)
        (directory / "root-report.json").write_text(json.dumps(report, indent=2) + "\n")
        if report["convex_metric_certified"]:
            return directory / "root-certificate.json", report
    except (AssertionError, ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        errors.append("Cartesian metric root: " + str(error))
    try:
        import hull_proposal

        hull, hull_report = hull_proposal.propose(polygons, numerical, p, digits)
        (directory / "hull-certificate.json").write_text(json.dumps(hull) + "\n")
        return directory / "hull-certificate.json", hull_report
    except (
        AssertionError,
        ValueError,
        ArithmeticError,
        ImportError,
        np.linalg.LinAlgError,
    ) as error:
        errors.append("Convex-hull metric: " + str(error))
    (directory / "convexity-limitations.json").write_text(json.dumps(errors, indent=2) + "\n")
    if report is not None:
        return directory / "root-certificate.json", report
    raise ValueError("; ".join(errors))


def produce(directory, digits):
    """Try the original metric, then the BI trace, without repeating proof replay.

    Route verification selects a proposal only. The launcher independently
    replays the selected certificate in a standard-library-only process.
    """
    directory = Path(directory)
    source = json.loads((directory / "input.json").read_text())
    numerical = json.loads((directory / "numerical.json").read_text())[0]
    original = dict(numerical, flips=[], triangles=m.initial_metric(source["faces"])[0])
    routes = [("original-metric", original)]
    if numerical.get("flips"):
        routes.append(("bi-trace", numerical))
    attempts, selected = [], None
    for name, estimate in routes:
        work = directory / ("proposal-" + name)
        work.mkdir()
        try:
            certificate, report = propose_metric(work, source["faces"], estimate, digits)
            attempts.append(dict(route=name, status=report["status"]))
            if selected is None or report["convex_metric_certified"]:
                selected = certificate
            if report["convex_metric_certified"]:
                break
        except (
            AssertionError,
            ValueError,
            ArithmeticError,
            ImportError,
            np.linalg.LinAlgError,
        ) as error:
            attempts.append(dict(route=name, error=str(error)))
    (directory / "proposal-routes.json").write_text(json.dumps(attempts, indent=2) + "\n")
    if selected is None:
        raise ValueError("No verified certificate proposal: " + json.dumps(attempts))
    shutil.copyfile(selected, directory / selected.name)
    (directory / "chosen-certificate.txt").write_text(selected.name + "\n")
    return selected.name


if __name__ == "__main__":
    import sys

    print(produce(Path(sys.argv[1]), int(sys.argv[2])))
