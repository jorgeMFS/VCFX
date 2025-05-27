#!/usr/bin/env bash

set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
BUILD_DIR="${ROOT_DIR}/build/python_bindings"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -DPYTHON_BINDINGS=ON ../..
make -j

cd "$SCRIPT_DIR"
source "$ROOT_DIR/add_vcfx_tools_to_path.sh"

PYTHONPATH="${BUILD_DIR}/python" python3 - <<'PY'
import vcfx
out = vcfx.trim("  hello  ")
if out != "hello":
    raise SystemExit('trim failed')
compressed = vcfx.read_maybe_compressed(b"hello")
if compressed != b"hello":
    raise SystemExit('read_maybe_compressed failed')
print('Python bindings OK:', out)
PY
