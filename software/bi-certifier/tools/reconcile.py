"""Compare completed fresh results with the saved corpus after all solving."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census-engine", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.census_engine.resolve()))
    from link_growth import canonical_code

    reference = {}
    for row in json.loads(args.corpus.read_text()):
        cert = json.loads(Path(row["path"]).read_text())
        code = canonical_code(cert["polygons"])
        if code in reference:
            raise ValueError("Duplicate reference complex")
        reference[code] = row["status"]
    rows = []
    seen = set()
    for path in sorted(args.results.glob("*/result.json")):
        result = json.loads(path.read_text())
        faces = json.loads((path.parent / "input.json").read_text())["faces"]
        code = canonical_code(faces)
        duplicate = code in seen
        seen.add(code)
        rows.append(
            dict(
                cid=result["cid"],
                status=result["status"],
                expected=reference.get(code),
                input_match=code in reference,
                duplicate=duplicate,
                verdict_match=result["status"] == reference.get(code),
                directory=str(path.parent),
            )
        )
    report = dict(
        count=len(rows),
        expected_count=len(reference),
        all_inputs_match=seen == set(reference) and not any(r["duplicate"] for r in rows),
        all_verdicts_match=all(r["verdict_match"] for r in rows),
        counts=dict(Counter(r["status"] for r in rows)),
        rows=rows,
    )
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2))
    return int(not (report["all_inputs_match"] and report["all_verdicts_match"]))


if __name__ == "__main__":
    sys.exit(main())
