# Verified corrections to numerical sphere guides

These files are inputs to the census launcher, not exact answers consumed by
the algebraic solver. Each correction is bound by SHA-256 to the original
partition, BI coordinate file, face file, and packaged replacement.

Version 3.0.9 includes one correction: census 66, CRP0-H97350c77009d16a.
Its old partition uniquely places the common center at vertex 7, producing an
exact zero-radius singleton. The replacement has three positive-radius groups
of sizes 9, 4, and 2, with a strictly interior center.

`run_body.py` applies this correction automatically when selecting the original
saved partition, including through `--partitions`. All other supplied partitions
retain their content. The original census files are not modified. The lower-level
`solve.py` and `run.sh` honor an explicitly supplied `--spheres` file as written.

The exact solver still derives its field and radii and verifies the original
body. The corresponding v3.0.9 results, independent check, and manual derivation
are distributed under `VALIDATION/v3.0.9-body66/`.
