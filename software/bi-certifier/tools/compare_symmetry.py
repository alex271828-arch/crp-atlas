"""Compare the direct map algorithm with v3's independent graph matcher."""

import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import symmetry

old_path = Path("/Users/ferleger/Downloads/crp_bi_numeric_v3/symmetry.py")
spec = importlib.util.spec_from_file_location("old_symmetry", old_path)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
rows = []
for item in json.loads((ROOT / "validation/corpus.json").read_text()):
    faces = json.loads(Path(item["path"]).read_text())["polygons"]
    start = time.monotonic()
    a = old.compute(faces)
    old_seconds = time.monotonic() - start
    start = time.monotonic()
    b = symmetry.compute(faces)
    new_seconds = time.monotonic() - start
    match = a["permutations"] == b["permutations"] and a["vertex_orbits"] == b["vertex_orbits"]
    rows.append(
        dict(
            case=Path(item["path"]).parent.name,
            match=match,
            order=b["order"],
            old_seconds=old_seconds,
            new_seconds=new_seconds,
        )
    )
    print(len(rows), match, round(old_seconds, 4), round(new_seconds, 4), flush=True)
(ROOT / "validation/symmetry-comparison.json").write_text(
    json.dumps(
        dict(
            all_match=all(r["match"] for r in rows),
            count=len(rows),
            old_seconds=sum(r["old_seconds"] for r in rows),
            new_seconds=sum(r["new_seconds"] for r in rows),
            rows=rows,
        ),
        indent=2,
    )
    + "\n"
)
