#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
VCFX_BIN="$ROOT_DIR/build/src/vcfx_wrapper/vcfx"

# Ensure built tools are in PATH so the wrapper can locate them
source "$ROOT_DIR/add_vcfx_tools_to_path.sh"

if [ ! -x "$VCFX_BIN" ]; then
  echo "vcfx executable not found: $VCFX_BIN"
  exit 1
fi

LIST_LONG="$($VCFX_BIN --list)"
LIST_ALIAS="$($VCFX_BIN list)"
if [ "$LIST_LONG" != "$LIST_ALIAS" ]; then
  echo "Output of 'vcfx list' does not match '--list'"
  diff <(echo "$LIST_LONG") <(echo "$LIST_ALIAS") || true
  exit 1
fi

echo "$LIST_LONG" > /dev/null # quiet shellcheck complaining about unused var

DOC_FIRST_LINE="$($VCFX_BIN help allele_counter | head -n 1)"
if ! echo "$DOC_FIRST_LINE" | grep -q "VCFX_allele_counter"; then
  echo "Help output for allele_counter does not show documentation"
  echo "First line was: $DOC_FIRST_LINE"
  exit 1
fi

echo "✓ vcfx wrapper tests passed"
