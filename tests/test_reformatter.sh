#!/usr/bin/env bash

set -e
set -o pipefail

###############################################################################
# Paths & Setup
###############################################################################
# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the root directory of the project
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Paths (work from either root or tests directory)
TOOL_BIN="${ROOT_DIR}/build/src/VCFX_reformatter/VCFX_reformatter"
TEST_DATA_DIR="${SCRIPT_DIR}/data/reformatter"
EXPECTED_DIR="${SCRIPT_DIR}/expected/reformatter"
OUTPUT_DIR="${SCRIPT_DIR}/out/reformatter"

# Create directories if they don't exist
mkdir -p "$TEST_DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$OUTPUT_DIR"

# Check if executable exists
if [ ! -f "$TOOL_BIN" ]; then
  echo "Error: $TOOL_BIN not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

###############################################################################
# Utility: convert literal "\t" to real tabs
###############################################################################
function makeTabs() {
    local inFile="$1"
    local outFile="$2"
    sed 's/\\t/\t/g' "$inFile" > "$outFile"
}

echo "=== Testing VCFX_reformatter ==="

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
HELP_OUT=$("$TOOL_BIN" --help 2>&1 || true)

if ! echo "$HELP_OUT" | grep -q "VCFX_reformatter: Reformat INFO/FORMAT fields"; then
    echo "✗ Test 1 failed: usage message not found."
    echo "Output was:"
    echo "$HELP_OUT"
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: No arguments => show help
###############################################################################
echo "Test 2: No arguments => usage"
ERR2=$("$TOOL_BIN" 2>&1 || true)

# The code calls displayHelp() if (argc==1)
if ! echo "$ERR2" | grep -q "Usage:"; then
    echo "✗ Test 2 failed: expected usage output for no arguments not found."
    echo "$ERR2"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Minimal VCF => compress-info with a single key
###############################################################################
echo "Test 3: Minimal VCF => compress-info"

# We'll create a minimal 2-line VCF with an INFO field including "AF=0.1;DP=10;FOO=xyz"
# We'll compress "DP" => remove it from INFO

cat > "$TEST_DATA_DIR/minimal_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\tPASS\tAF=0.1;DP=10;FOO=xyz
1\t200\t.\tT\tC\t.\tPASS\tDP=20;AF=0.5
EOF

makeTabs "$TEST_DATA_DIR/minimal_tmp.vcf" "$TEST_DATA_DIR/minimal.vcf"

OUT3=$("$TOOL_BIN" --compress-info "DP" < "$TEST_DATA_DIR/minimal.vcf")

# Save the output for reference
echo "$OUT3" > "$OUTPUT_DIR/compressed_output.vcf"

# Check that the lines with #CHROM are present
if ! echo "$OUT3" | grep -q "^#CHROM"; then
    echo "✗ Test 3 failed: #CHROM header line missing."
    echo "$OUT3"
    exit 1
fi

# We expect:
#   line at POS=100 => "AF=0.1;FOO=xyz" (DP removed)
#   line at POS=200 => "AF=0.5" (DP removed)
if ! echo "$OUT3" | grep -qE "^1\s+100.*AF=0\.1;FOO=xyz"; then
    echo "✗ Test 3 failed: expected 'AF=0.1;FOO=xyz' in POS=100 line."
    echo "$OUT3"
    exit 1
fi
if ! echo "$OUT3" | grep -qE "^1\s+200.*AF=0\.5" ; then
    echo "✗ Test 3 failed: expected 'AF=0.5' in POS=200 line after removing DP."
    echo "$OUT3"
    exit 1
fi

echo "✓ Test 3 passed"

###############################################################################
# Test 4: Minimal VCF => reorder-info (like "FOO,AF" => FOO=...,AF=..., leftover appended)
###############################################################################
echo "Test 4: reorder-info"

cat > "$TEST_DATA_DIR/reorder_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\tPASS\tAF=0.1;DP=10;FOO=xyz
EOF

makeTabs "$TEST_DATA_DIR/reorder_tmp.vcf" "$TEST_DATA_DIR/reorder.vcf"

# We'll reorder so that FOO is first, then AF, leftover => DP
OUT4=$("$TOOL_BIN" --reorder-info "FOO,AF" < "$TEST_DATA_DIR/reorder.vcf")

# Save the output for reference
echo "$OUT4" > "$OUTPUT_DIR/reordered_output.vcf"

# Expect => "FOO=xyz;AF=0.1;DP=10"
if ! echo "$OUT4" | grep -q "FOO=xyz;AF=0.1;DP=10"; then
    echo "✗ Test 4 failed: reorder-info did not produce 'FOO=xyz;AF=0.1;DP=10'."
    echo "$OUT4"
    exit 1
fi

echo "✓ Test 4 passed"

###############################################################################
# Test 5: line with <8 columns => skipping
###############################################################################
echo "Test 5: skip lines with <8 columns"

cat > "$TEST_DATA_DIR/fewcols_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\t.\tPASS\t
EOF

makeTabs "$TEST_DATA_DIR/fewcols_tmp.vcf" "$TEST_DATA_DIR/fewcols.vcf"

# The single data line has only 7 columns => warning => skip
OUT5=$("$TOOL_BIN" --compress-info "DP" < "$TEST_DATA_DIR/fewcols.vcf" 2>&1 || true)

# Save the output for reference
echo "$OUT5" > "$OUTPUT_DIR/fewcols_output.txt"

# Check that we see a "Warning: line with <8 columns => skipping."
if ! echo "$OUT5" | grep -q "Warning: line with <8 columns => skipping."; then
    echo "✗ Test 5 failed: expected warning about <8 columns not found."
    echo "$OUT5"
    exit 1
fi

# The output should still have the #CHROM line
if ! echo "$OUT5" | grep -q "^#CHROM"; then
    echo "✗ Test 5 failed: missing #CHROM line in output."
    echo "$OUT5"
    exit 1
fi

echo "✓ Test 5 passed"

###############################################################################
echo "✅ All tests for VCFX_reformatter passed successfully!"
