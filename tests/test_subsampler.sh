#!/usr/bin/env bash

set -e
set -o pipefail

###############################################################################
# Paths & Setup
###############################################################################
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

TOOL_BIN="${ROOT_DIR}/build/src/VCFX_subsampler/VCFX_subsampler"
TEST_DATA_DIR="${SCRIPT_DIR}/data/subsampler"
EXPECTED_DIR="${SCRIPT_DIR}/expected/subsampler"
OUTPUT_DIR="${SCRIPT_DIR}/out/subsampler"

mkdir -p "$TEST_DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$OUTPUT_DIR"

# Utility function: Convert literal "\t" into real tab characters
function makeTabs() {
    local inFile="$1"
    local outFile="$2"
    sed 's/\\t/\t/g' "$inFile" > "$outFile"
}

if [ ! -x "$TOOL_BIN" ]; then
  echo "Error: $TOOL_BIN not found or not executable."
  exit 1
fi

echo "=== Testing VCFX_subsampler ==="

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
HELP_OUT=$("$TOOL_BIN" --help 2>&1 || true)
if ! echo "$HELP_OUT" | grep -q "VCFX_subsampler: Randomly pick N lines from a VCF data section."; then
    echo "✗ Test 1 failed: usage message not found."
    echo "Output was:"
    echo "$HELP_OUT"
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: No arguments => usage
###############################################################################
echo "Test 2: No arguments => usage"
NOARGS_OUT=$("$TOOL_BIN" 2>&1 || true)
if ! echo "$NOARGS_OUT" | grep -q "VCFX_subsampler: Randomly pick N lines from a VCF"; then
    echo "✗ Test 2 failed: expected usage text when no arguments are given."
    echo "$NOARGS_OUT"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Invalid --subsample (negative) => error
###############################################################################
echo "Test 3: Invalid --subsample"
ERR3=$("$TOOL_BIN" --subsample "-5" 2>&1 || true)
if ! echo "$ERR3" | grep -q "Error: invalid subsample size."; then
    echo "✗ Test 3 failed: expected error about invalid subsample size."
    echo "$ERR3"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Minimal VCF => keep all valid lines if fewer than N, warn on invalid lines
###############################################################################
echo "Test 4: Minimal VCF => keep them all if fewer than N"

# Create a minimal VCF file using a single-quoted heredoc so that "\t" remains literal.
cat > "$TEST_DATA_DIR/minimal_tmp.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tSomething
chr1\tThisIsInvalid
chr1\t200\t.\tT\tC\t.\tPASS\tAnother
EOF

# Convert literal "\t" into actual tabs.
makeTabs "$TEST_DATA_DIR/minimal_tmp.vcf" "$TEST_DATA_DIR/minimal.vcf"

# Run subsampler with a sample size larger than the number of valid data lines.
# Capture both stdout and stderr.
OUT4=$("$TOOL_BIN" --subsample 5 < "$TEST_DATA_DIR/minimal.vcf" 2>&1 || true)
echo "$OUT4" > "$OUTPUT_DIR/minimal_output.txt"

# 1) The header line must appear.
if ! echo "$OUT4" | grep -q "^#CHROM"; then
    echo "✗ Test 4 failed: missing #CHROM header line."
    echo "$OUT4"
    exit 1
fi

# 2) Expect a warning about skipping a line with fewer than 8 columns.
if ! echo "$OUT4" | grep -q "Warning: skipping line with <8 columns."; then
    echo "✗ Test 4 failed: expected warning about <8 columns not found."
    echo "$OUT4"
    exit 1
fi

# 3) Expect the two valid data lines (POS=100 and POS=200) to appear.
if ! echo "$OUT4" | grep -q "^chr1[[:space:]]100"; then
    echo "✗ Test 4 failed: missing data line with pos=100."
    echo "$OUT4"
    exit 1
fi
if ! echo "$OUT4" | grep -q "^chr1[[:space:]]200"; then
    echo "✗ Test 4 failed: missing data line with pos=200."
    echo "$OUT4"
    exit 1
fi
echo "✓ Test 4 passed"

###############################################################################
# Test 5: Reservoir sampling with fixed seed => deterministic subset
###############################################################################
echo "Test 5: Reservoir sampling with fixed seed"

cat > "$TEST_DATA_DIR/seedsample_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tX
chr1\t200\t.\tT\tC\t.\tPASS\tX
chr1\t300\t.\tG\tA\t.\tPASS\tX
chr1\t400\t.\tC\tT\t.\tPASS\tX
EOF

makeTabs "$TEST_DATA_DIR/seedsample_tmp.vcf" "$TEST_DATA_DIR/seedsample.vcf"

# Run with --subsample 2 and a fixed seed (e.g. 42)
OUT5=$("$TOOL_BIN" --subsample 2 --seed 42 < "$TEST_DATA_DIR/seedsample.vcf")
echo "$OUT5" > "$OUTPUT_DIR/seedsample_output.txt"

# Verify header is present.
if ! echo "$OUT5" | grep -q "^#CHROM"; then
    echo "✗ Test 5 failed: missing #CHROM in output."
    echo "$OUT5"
    exit 1
fi

# Count the valid data lines (those not starting with '#')
dataCount=$(echo "$OUT5" | grep -E "^[^#]" | wc -l | tr -d ' ')
if [ "$dataCount" -ne 2 ]; then
    echo "✗ Test 5 failed: expected exactly 2 data lines, found $dataCount."
    echo "$OUT5"
    exit 1
fi

echo "✓ Test 5 passed"

###############################################################################
echo "✅ All tests for VCFX_subsampler passed successfully!"
