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
TOOL_BIN="${ROOT_DIR}/build/src/VCFX_metadata_summarizer/VCFX_metadata_summarizer"
TEST_DATA_DIR="${SCRIPT_DIR}/data/metadata_summarizer"
EXPECTED_DIR="${SCRIPT_DIR}/expected/metadata_summarizer"
OUTPUT_DIR="${SCRIPT_DIR}/out/metadata_summarizer"

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

echo "=== Testing VCFX_metadata_summarizer ==="

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
HELP_OUT=$("$TOOL_BIN" --help 2>&1 || true)

# Check if the help text mentions "VCFX_metadata_summarizer: Summarize key metadata"
if ! echo "$HELP_OUT" | grep -q "VCFX_metadata_summarizer: Summarize key metadata"; then
    echo "✗ Test 1 failed: help message not found."
    echo "Output:"
    echo "$HELP_OUT"
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: Minimal VCF => Check summary
###############################################################################
echo "Test 2: Minimal VCF"

# Create a minimal VCF with one contig, one INFO, and a single variant
cat > "$TEST_DATA_DIR/minimal.vcf" <<EOF
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tDP=10
EOF

# Run the summarizer
OUT2=$("$TOOL_BIN" < "$TEST_DATA_DIR/minimal.vcf")

# Save output for reference
echo "$OUT2" > "$OUTPUT_DIR/minimal_summary.txt"

# Check the lines for correctness:
# We expect:
#   "Number of unique contigs: 1"
#   "Number of unique INFO fields: 1"
#   "Number of unique FILTER fields: 0"
#   "Number of unique FORMAT fields: 0"
#   "Number of samples: 0"        (no sample columns beyond #CHROM..INFO)
#   "Number of variants: 1"

if ! echo "$OUT2" | grep -q "Number of unique contigs: 1"; then
    echo "✗ Test 2 failed: contigs count mismatch"
    echo "$OUT2"
    exit 1
fi

if ! echo "$OUT2" | grep -q "Number of unique INFO fields: 1"; then
    echo "✗ Test 2 failed: info fields count mismatch"
    echo "$OUT2"
    exit 1
fi

if ! echo "$OUT2" | grep -q "Number of samples: 0"; then
    echo "✗ Test 2 failed: samples mismatch"
    echo "$OUT2"
    exit 1
fi

if ! echo "$OUT2" | grep -q "Number of variants: 1"; then
    echo "✗ Test 2 failed: variants mismatch"
    echo "$OUT2"
    exit 1
fi

echo "✓ Test 2 passed"

###############################################################################
echo "✅ All tests for VCFX_metadata_summarizer passed successfully!"
