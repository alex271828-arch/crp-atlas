"""Plain-text coordinate export from an independently verified certificate."""

import json
from pathlib import Path


def render(folder):
    folder = Path(folder)
    verdict = json.loads((folder / "verification.json").read_text())
    if verdict.get("verified") is not True:
        raise ValueError("Cannot label unverified coordinates as certified")
    c = json.loads((folder / "coordinates.json").read_text())
    s = c["source_algebra"]
    lines = [
        "CERTIFIED EXACT VERTEX COORDINATES",
        "",
        f"Primitive symbol: {s['primitive_symbol']}",
        f"Defining polynomial: {s['primitive_polynomial']} = 0",
        f"Real isolating interval: {s['isolating_interval']}",
        "",
        "Positive square roots:",
    ]
    lines += [f"r{i} = positive sqrt({r})" for i, r in enumerate(c["positive_radicals"])]
    lines += ["", "Metric G:"] + ["  " + str(row) for row in c["metric"]]
    for title, key in [
        ("Metric coordinates (squared distance = delta^T G delta)", "metric_coordinate_formulas"),
        ("Cartesian coordinates (ordinary Euclidean metric)", "cartesian_coordinate_formulas"),
    ]:
        lines += ["", title] + [f"{i}: (" + ", ".join(v) + ")" for i, v in enumerate(c[key])]
    lines += [
        "",
        "Both representations are exactly verified. No minimal field or globally shortest expression is claimed.",
        "The full centers, isometry and original labeled faces are in coordinates.json.",
    ]
    with (folder / "coordinates.txt").open("x") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    import sys

    render(sys.argv[1])
