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
TOOL_BIN="${ROOT_DIR}/build/src/VCFX_phase_quality_filter/VCFX_phase_quality_filter"
TEST_DATA_DIR="${SCRIPT_DIR}/data/phase_quality_filter"
EXPECTED_DIR="${SCRIPT_DIR}/expected/phase_quality_filter"
OUTPUT_DIR="${SCRIPT_DIR}/out/phase_quality_filter"

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
# Utility function to convert literal "\t" into actual tab characters
###############################################################################
function makeTabs() {
    local inFile="$1"
    local outFile="$2"
    sed 's/\\t/\t/g' "$inFile" > "$outFile"
}

echo "=== Testing VCFX_phase_quality_filter ==="

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
HELP_OUT=$("$TOOL_BIN" --help 2>&1 || true)

if ! echo "$HELP_OUT" | grep -q "VCFX_phase_quality_filter: Filter variants by phasing quality"; then
    echo "✗ Test 1 failed: help message not found."
    echo "Output was:"
    echo "$HELP_OUT"
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: Missing --filter-pq => expect an error
###############################################################################
echo "Test 2: Missing --filter-pq"
ERR2=$("$TOOL_BIN" 2>&1 || true)

if ! echo "$ERR2" | grep -q "Error: Must specify condition with --filter-pq"; then
    echo "✗ Test 2 failed: expected error about missing --filter-pq not found."
    echo "$ERR2"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Minimal VCF => check which lines are kept by PQ>30
###############################################################################
echo "Test 3: Minimal VCF => test keep/drops"

cat > "$TEST_DATA_DIR/minimal_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\tPASS\tPQ=5
1\t200\t.\tT\tC\t.\tPASS\tPQ=30
1\t300\t.\tG\tA\t.\tPASS\tPQ=31
1\t400\t.\tC\tT\t.\tPASS\t.   # Missing => PQ=0
1\t500\t.\tA\tG\t.\tPASS\tPQ=20
EOF

# Convert literal "\t" => real tabs
makeTabs "$TEST_DATA_DIR/minimal_tmp.vcf" "$TEST_DATA_DIR/minimal.vcf"

# Now run the filter with PQ>30 => we should keep only line(s) that have PQ>30
OUT3=$("$TOOL_BIN" --filter-pq "PQ>30" < "$TEST_DATA_DIR/minimal.vcf")

# Save the output for reference
echo "$OUT3" > "$OUTPUT_DIR/filtered_output.vcf"

# We expect just the #CHROM line plus the data line with POS=300
# 1) Check that POS=300 is present
if ! echo "$OUT3" | grep -q "^1[[:space:]]300"; then
    echo "✗ Test 3 failed: expected only POS=300 line to remain."
    echo "$OUT3"
    exit 1
fi

# 2) POS=200 (PQ=30) should be filtered out
if echo "$OUT3" | grep -q "^1[[:space:]]200"; then
    echo "✗ Test 3 failed: line with PQ=30 should be filtered out."
    exit 1
fi

# 3) Check for the #CHROM header line
if ! echo "$OUT3" | grep -q "^#CHROM"; then
    echo "✗ Test 3 failed: missing #CHROM line in output."
    echo "$OUT3"
    exit 1
fi

echo "✓ Test 3 passed"

###############################################################################
# Test 4: Invalid condition => parse error
###############################################################################
echo "Test 4: Invalid condition"

ERR4=$("$TOOL_BIN" --filter-pq "PQ@100" < /dev/null 2>&1 || true)

if ! echo "$ERR4" | grep -q "Unable to parse condition 'PQ@100'"; then
    echo "✗ Test 4 failed: expected parse error not found."
    echo "$ERR4"
    exit 1
fi
echo "✓ Test 4 passed"

###############################################################################
echo "✅ All tests for VCFX_phase_quality_filter passed successfully!"
