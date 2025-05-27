#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

VCFX_BUILD_BIN="$ROOT_DIR/build/src/vcfx_wrapper/vcfx"
TMPDIR="$(mktemp -d)"
mkdir -p "$TMPDIR/bin"
cp "$VCFX_BUILD_BIN" "$TMPDIR/bin/vcfx"
# copy docs to share/doc/VCFX relative to bin
mkdir -p "$TMPDIR/share/doc/VCFX"
cp -r "$ROOT_DIR/docs/." "$TMPDIR/share/doc/VCFX"

PATH="$TMPDIR/bin:$PATH" vcfx help allele_counter > "$TMPDIR/out.txt"
FIRST_LINE="$(head -n 1 "$TMPDIR/out.txt")"
rm -rf "$TMPDIR"

if ! echo "$FIRST_LINE" | grep -q "VCFX_allele_counter"; then
  echo "Documentation lookup failed"
  exit 1
fi

echo "✓ doc lookup portable test passed"
