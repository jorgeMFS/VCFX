#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

TMPDIR="$(mktemp -d)"

# create dummy executables
cat <<'SH' > "$TMPDIR/VCFX_dummy1"
#!/bin/sh
exit 0
SH
chmod +x "$TMPDIR/VCFX_dummy1"

cat <<'SH' > "$TMPDIR/VCFX_dummy2"
#!/bin/sh
exit 0
SH
chmod +x "$TMPDIR/VCFX_dummy2"

# limit PATH to ensure vcfx is not found, but python remains available
export PATH="$TMPDIR:/usr/bin:/bin"

PYTHONPATH="$ROOT_DIR/build/python" python3 - <<'PY'
import vcfx
expected = {"dummy1", "dummy2"}
result = set(vcfx.available_tools())
if result != expected:
    raise SystemExit(f"unexpected tools: {result}")
print("✓ available_tools fallback passed")
PY

rm -rf "$TMPDIR"

