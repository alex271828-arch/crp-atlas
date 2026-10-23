"""Join algebraic source cases to the BI 4.0 survivors by full canonical incidence."""

from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, "/Users/ferleger/Downloads/crp3456810_census_v1_16")
from link_growth import canonical_code

batch = Path(sys.argv[1])
manifest = json.loads((batch / "manifest.json").read_text())
bi = Path("/Users/ferleger/Downloads/crp_bi_certifier_v4/validation/final-census")
reference = {}
for p in bi.glob("*/result.json"):
    d = json.loads(p.read_text())
    if d["status"] != "OPEN_FOR_ALGEBRAIC":
        continue
    faces = json.loads((p.parent / "input.json").read_text())["faces"]
    code = canonical_code(faces)
    if code in reference:
        raise ValueError("Duplicate BI survivor")
    reference[code] = dict(bi_case=p.parent.name[:3], bi_cid=d["cid"], bi_directory=str(p.parent))
rows = []
seen = set()
for item in manifest["sources"]:
    src = Path(item["source"])
    old = Path(item["original_batch_source"])
    case = json.loads((src / "input.json").read_text())
    code = canonical_code(case["faces"])
    if code not in reference or code in seen:
        raise ValueError("Missing or duplicate input complex")
    original = json.loads((old / "input.json").read_text())
    if canonical_code(original["faces"]) != code:
        raise ValueError("Repair changed original face complex")
    seen.add(code)
    rows.append(
        dict(
            **item,
            **reference[code],
            faces_match=True,
            repair_incidence_matches=True,
            source_hashes={
                n: hashlib.sha256((src / n).read_bytes()).hexdigest()
                for n in ["input.json", "radii.json"]
            }
        )
    )
assert len(rows) == 198 and seen == set(reference)
report = dict(count=198, all_bi_survivors_matched=True, all_repairs_match=True, rows=rows)
with (batch / "incidence-ledger.json").open("x") as f:
    json.dump(report, f, indent=2)
    f.write("\n")
print("Matched all 198 source complexes, including repairs, to all 198 BI 4.0 survivors.")
