"""Convert the 195 saved successes plus the three matched repairs, in isolation."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import json, subprocess, sys, time, hashlib
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
source = Path("/Users/ferleger/Downloads/crp_algebraic_radii_v3/out")
batch = source / "all-198-v3.0.7-20260916-050904"
output = Path(sys.argv[1]).resolve()
output.mkdir(parents=True, exist_ok=False)
items = []
for p in sorted(batch.glob("*/result.json")):
    d = json.loads(p.read_text())
    number = p.parent.name[:3]
    original = p.parent
    if d["status"] != "SOLVED":
        p = next(source.glob(number + "-CRP*/result.json"))
        d = json.loads(p.read_text())
    if not d["certified"]:
        raise ValueError("Uncertified source")
    items.append(
        dict(
            case=number,
            source=str(p.parent),
            original_batch_source=str(original),
            version=d["version"],
        )
    )
assert len(items) == 198
manifest = dict(
    sources=items,
    code_sha256={
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [*ROOT.glob("*.py"), *ROOT.glob("vendor/*.py")]
    },
)
(output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def one(item):
    dest = output / item["case"]
    start = time.monotonic()
    p = subprocess.run(
        [
            sys.executable,
            "-B",
            str(ROOT / "convert.py"),
            item["source"],
            "--output",
            str(dest),
            "--timeout",
            "300",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    (output / (item["case"] + ".log")).write_text(p.stdout)
    d = json.loads((dest / "result.json").read_text())
    return dict(**item, **d, wall_seconds=time.monotonic() - start, returncode=p.returncode)


rows = []
start = time.monotonic()
with ThreadPoolExecutor(max_workers=3) as pool:
    for future in as_completed([pool.submit(one, item) for item in items]):
        row = future.result()
        rows.append(row)
        print(len(rows), row["case"], row["status"], round(row["seconds"], 2), flush=True)
        (output / "summary.json").write_text(
            json.dumps(sorted(rows, key=lambda x: x["case"]), indent=2) + "\n"
        )
unchanged = all(
    hashlib.sha256((ROOT / n).read_bytes()).hexdigest() == h
    for n, h in manifest["code_sha256"].items()
)
(output / "completion.json").write_text(
    json.dumps(
        dict(
            total=len(rows),
            counts=dict(Counter(r["status"] for r in rows)),
            seconds=time.monotonic() - start,
            source_unchanged=unchanged,
        ),
        indent=2,
    )
    + "\n"
)
