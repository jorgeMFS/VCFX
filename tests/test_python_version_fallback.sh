#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Determine package version from pyproject.toml
VERSION=$(grep '^version' "$ROOT_DIR/python/pyproject.toml" | cut -d '"' -f2)
TMPDIR="$(mktemp -d)"
mkdir -p "$TMPDIR/vcfx"
cp "$ROOT_DIR"/python/*.py "$TMPDIR/vcfx"
cp "$ROOT_DIR"/python/py.typed "$TMPDIR/vcfx"
mkdir -p "$TMPDIR/vcfx-$VERSION.dist-info"
cat > "$TMPDIR/vcfx-$VERSION.dist-info/METADATA" <<EOF2
Metadata-Version: 2.1
Name: vcfx
Version: $VERSION
EOF2

# Provide a dummy vcfx wrapper to satisfy available_tools()
mkdir -p "$TMPDIR/bin"
cat <<'SH' > "$TMPDIR/bin/vcfx"
#!/bin/sh
exit 0
SH
chmod +x "$TMPDIR/bin/vcfx"
export PATH="$TMPDIR/bin:/usr/bin:/bin"

PYTHONPATH="$TMPDIR" python3 - <<PY
import vcfx
expected = "${VERSION}"
if vcfx.__version__ != expected:
    raise SystemExit(f"expected {expected}, got {vcfx.__version__}")
print("\u2713 version fallback passed")
PY

rm -rf "$TMPDIR"
