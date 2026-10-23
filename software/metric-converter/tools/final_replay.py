"""Verify the complete coordinate collection with the frozen release checker."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from collections import Counter
import hashlib, json, shutil, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=False)
ledger = json.loads((source / "incidence-ledger.json").read_text())
rows = ledger["rows"]
costs = {r["case"]: r["seconds"] for r in json.loads((source / "summary.json").read_text())}
rows = sorted(rows, key=lambda r: costs.get(r["case"], 0), reverse=True)
assert len(rows) == 198
hashes = {
    str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in [*ROOT.glob("*.py"), *ROOT.glob("vendor/*.py")]
}
(out / "manifest.json").write_text(
    json.dumps(
        dict(code_sha256=hashes, source_ledger=str(source / "incidence-ledger.json"), sources=rows),
        indent=2,
    )
    + "\n"
)
for row in rows:
    original = Path(row["source"])
    cert = json.loads((source / row["case"] / "coordinates.json").read_text())
    case = json.loads((original / "input.json").read_text())
    saved = json.loads((original / "radii.json").read_text())
    assert cert["input"] == {k: case[k] for k in ["faces", "groups"]}
    assert cert["source_algebra"] == {
        k: saved[k]
        for k in [
            "primitive_symbol",
            "primitive_polynomial",
            "isolating_interval",
            "groups",
            "squared_radii",
        ]
    }
    for name, h in row["source_hashes"].items():
        assert (
            hashlib.sha256((original / name).read_bytes()).hexdigest()
            == h
            == cert["source_files"][name]["sha256"]
        )
    target = out / row["case"]
    target.mkdir()
    shutil.copy2(source / row["case"] / "coordinates.json", target / "coordinates.json")
start = time.monotonic()


def one(row):
    target = out / row["case"]
    t = time.monotonic()
    with (target / "check.log").open("x") as log:
        p = subprocess.run(
            [
                sys.executable,
                "-B",
                str(ROOT / "check.py"),
                str(target / "coordinates.json"),
                "--output",
                str(target / "verification.json"),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=600,
        )
    r = json.loads((target / "verification.json").read_text())
    assert p.returncode == 0 and r["verified"], (row["case"], r)
    sys.path.insert(0, str(ROOT))
    from report import render

    render(target)
    return dict(
        case=row["case"],
        bi_case=row["bi_case"],
        bi_cid=row["bi_cid"],
        source_version=row["version"],
        source=row["source"],
        certificate_sha256=hashlib.sha256((target / "coordinates.json").read_bytes()).hexdigest(),
        seconds=time.monotonic() - t,
        **r
    )


results = []
with ThreadPoolExecutor(max_workers=6) as pool:
    for f in as_completed([pool.submit(one, row) for row in rows]):
        r = f.result()
        results.append(r)
        print(len(results), r["case"], r["status"], round(r["seconds"], 2), flush=True)
        (out / "summary.json").write_text(
            json.dumps(sorted(results, key=lambda r: r["case"]), indent=2) + "\n"
        )
unchanged = all(hashlib.sha256((ROOT / n).read_bytes()).hexdigest() == h for n, h in hashes.items())
assert unchanged
(out / "completion.json").write_text(
    json.dumps(
        dict(
            total=len(results),
            counts=dict(Counter(r["status"] for r in results)),
            source_unchanged=unchanged,
            all_inputs_matched=True,
            seconds=time.monotonic() - start,
        ),
        indent=2,
    )
    + "\n"
)
