#!/bin/sh
set -eu
crp_script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ -n "${CRP_RADII_PYTHON:-}" ]; then
    crp_python=$CRP_RADII_PYTHON
elif [ -x /Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10 ]; then
    crp_python=/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10
else
    crp_python=python3
fi
exec "$crp_python" "$crp_script_dir/solve.py" "$@"
