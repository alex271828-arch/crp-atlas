#!/bin/sh
set -eu
cd "$(dirname "$0")"
# Use the separately installed Command Line Tools if the selected Xcode cannot
# run. This changes only this build's environment, not system settings/licenses.
if [ "$(uname -s)" = Darwin ] && [ -z "${DEVELOPER_DIR:-}" ] &&
   ! /usr/bin/xcrun --find clang++ >/dev/null 2>&1 &&
   [ -x /Library/Developer/CommandLineTools/usr/bin/clang++ ]; then
  DEVELOPER_DIR=/Library/Developer/CommandLineTools
  export DEVELOPER_DIR
fi
CMAKE_BIN="${CMAKE:-cmake}"
if ! command -v "$CMAKE_BIN" >/dev/null 2>&1; then
  sage_cmake=/Applications/SageMath-10-0.app/Contents/Frameworks/Sage.framework/Versions/10.0/local/bin/cmake
  if [ -x "$sage_cmake" ]; then CMAKE_BIN="$sage_cmake"; else
    echo "CMake is required; install it or set CMAKE to its executable." >&2
    exit 1
  fi
fi
"$CMAKE_BIN" -S . -B build -DCMAKE_BUILD_TYPE=Release -DSUITESPARSE=OFF
"$CMAKE_BIN" --build build --target metric_probe gradient_check -j "${BUILD_JOBS:-2}"
cp build/metric_probe ./metric_probe
cp build/gradient_check ./gradient_check
