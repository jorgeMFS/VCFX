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
TOOL_BIN="${ROOT_DIR}/build/src/VCFX_region_subsampler/VCFX_region_subsampler"
TEST_DATA_DIR="${SCRIPT_DIR}/data/region_subsampler"
EXPECTED_DIR="${SCRIPT_DIR}/expected/region_subsampler"
OUTPUT_DIR="${SCRIPT_DIR}/out/region_subsampler"

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
# Utility: Convert literal "\t" to real ASCII tabs
###############################################################################
function makeTabs() {
    local inFile="$1"
    local outFile="$2"
    sed 's/\\t/\t/g' "$inFile" > "$outFile"
}

echo "=== Testing VCFX_region_subsampler ==="

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
HELP_OUT=$("$TOOL_BIN" --help 2>&1 || true)

# We expect a line like "VCFX_region_subsampler: Keep only variants ..."
if ! echo "$HELP_OUT" | grep -q "VCFX_region_subsampler: Keep only variants whose (CHROM,POS)"; then
    echo "✗ Test 1 failed: help message not found."
    echo "Output was:"
    echo "$HELP_OUT"
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: No arguments => expect usage error about region-bed
###############################################################################
echo "Test 2: No arguments => usage error"
ERR2=$("$TOOL_BIN" 2>&1 || true)

# The code prints "Error: Must specify --region-bed <FILE>."
if ! echo "$ERR2" | grep -q "Error: Must specify --region-bed <FILE>."; then
    echo "✗ Test 2 failed: expected error about missing --region-bed not found."
    echo "$ERR2"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Minimal BED => keep only [chr1:1-100]
###############################################################################
echo "Test 3: Minimal BED => keep region [chr1:1..100]"

# Create a minimal BED with a single region => 0-based start=0, end=100 => final region is [1..100]
cat > "$TEST_DATA_DIR/regionA.bed" <<EOF
chr1	0	100
EOF

# Then create a minimal VCF with variants at pos=50, pos=101 => only pos=50 should remain
cat > "$TEST_DATA_DIR/minimal_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t50\t.\tA\tG\t.\tPASS\tSOMETHING
chr1\t101\t.\tT\tC\t.\tPASS\tANOTHER
chr2\t50\t.\tG\tA\t.\tPASS\tOUTSIDE
EOF

makeTabs "$TEST_DATA_DIR/minimal_tmp.vcf" "$TEST_DATA_DIR/minimal.vcf"

OUT3=$("$TOOL_BIN" --region-bed "$TEST_DATA_DIR/regionA.bed" < "$TEST_DATA_DIR/minimal.vcf")

# Save the output for reference
echo "$OUT3" > "$OUTPUT_DIR/minimal_output.vcf"

# We expect:
#  - #CHROM line present
#  - line with pos=50 remains
#  - line with pos=101 is gone
#  - line with chr2 => gone

# 1) Check if #CHROM is present
if ! echo "$OUT3" | grep -q "^#CHROM"; then
    echo "✗ Test 3 failed: missing header line (#CHROM)."
    echo "$OUT3"
    exit 1
fi

# 2) check for pos=50 => must remain
if ! echo "$OUT3" | grep -q "^chr1[[:space:]]50"; then
    echo "✗ Test 3 failed: line with chr1 pos=50 not found."
    echo "$OUT3"
    exit 1
fi

# 3) check pos=101 => must be absent
if echo "$OUT3" | grep -q "^chr1[[:space:]]101"; then
    echo "✗ Test 3 failed: line with pos=101 should be filtered out."
    exit 1
fi

# 4) check chr2 => absent
if echo "$OUT3" | grep -q "^chr2"; then
    echo "✗ Test 3 failed: line with chr2 pos=50 is outside region => should be filtered."
    exit 1
fi

echo "✓ Test 3 passed"

###############################################################################
# Test 4: Multi-line BED => keep lines that fall in any region
###############################################################################
echo "Test 4: Multi-line BED => multiple intervals"

cat > "$TEST_DATA_DIR/regionB.bed" <<EOF
chr1	0	100
chr2	100	200
EOF

# This means we keep:
#   - [chr1:1..100]
#   - [chr2:101..200]
cat > "$TEST_DATA_DIR/multi_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tIN_RANGE1
chr1\t101\t.\tT\tC\t.\tPASS\tOUT_RANGE1
chr2\t101\t.\tG\tA\t.\tPASS\tIN_RANGE2
chr2\t200\t.\tC\tT\t.\tPASS\tIN_RANGE2B
chr2\t201\t.\tA\tC\t.\tPASS\tOUT_RANGE2
EOF

makeTabs "$TEST_DATA_DIR/multi_tmp.vcf" "$TEST_DATA_DIR/multi.vcf"

OUT4=$("$TOOL_BIN" --region-bed "$TEST_DATA_DIR/regionB.bed" < "$TEST_DATA_DIR/multi.vcf")

# Save the output for reference
echo "$OUT4" > "$OUTPUT_DIR/multi_output.vcf"

# Expect:
#   - keep chr1 pos=100
#   - skip chr1 pos=101
#   - keep chr2 pos=101, pos=200
#   - skip chr2 pos=201
#   - keep #CHROM

# 1) must have #CHROM
if ! echo "$OUT4" | grep -q "^#CHROM"; then
    echo "✗ Test 4 failed: missing #CHROM line."
    echo "$OUT4"
    exit 1
fi

# 2) check each data line
if ! echo "$OUT4" | grep -q "^chr1[[:space:]]100"; then
    echo "✗ Test 4 failed: expected chr1:100 to be kept."
    echo "$OUT4"
    exit 1
fi
if echo "$OUT4" | grep -q "^chr1[[:space:]]101"; then
    echo "✗ Test 4 failed: chr1:101 should be out of region => skip."
    exit 1
fi
if ! echo "$OUT4" | grep -q "^chr2[[:space:]]101"; then
    echo "✗ Test 4 failed: expected chr2:101 to remain."
    echo "$OUT4"
    exit 1
fi
if ! echo "$OUT4" | grep -q "^chr2[[:space:]]200"; then
    echo "✗ Test 4 failed: expected chr2:200 to remain."
    echo "$OUT4"
    exit 1
fi
if echo "$OUT4" | grep -q "^chr2[[:space:]]201"; then
    echo "✗ Test 4 failed: chr2:201 out of region => skip."
    exit 1
fi

echo "✓ Test 4 passed"

###############################################################################
# Test 5: BED file with invalid lines => skip them
###############################################################################
echo "Test 5: BED file with invalid lines"

# 1) Create a bad BED with an invalid line and a negative interval
cat > "$TEST_DATA_DIR/bad_bed.bed" <<EOF
chr1	0	100
thisIsInvalid
chr2	50	40   # start> end => skip
EOF

# 2) Create a mini VCF (with literal \t) that has one line that should pass
cat > "$TEST_DATA_DIR/bad_bed_vcf_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tOK
EOF

# Convert literal "\t" => real tabs
makeTabs "$TEST_DATA_DIR/bad_bed_vcf_tmp.vcf" "$TEST_DATA_DIR/bad_bed_vcf.vcf"

# 3) Run the tool: region=chr1:1..100 => pos=100 is in range, so we expect it to remain
OUT5=$("$TOOL_BIN" --region-bed "$TEST_DATA_DIR/bad_bed.bed" < "$TEST_DATA_DIR/bad_bed_vcf.vcf" 2>&1 || true)

# Save the output for reference
echo "$OUT5" > "$OUTPUT_DIR/bad_bed_output.txt"

# The tool should warn about the invalid bed line "thisIsInvalid"
if ! echo "$OUT5" | grep -q "Warning: skipping invalid bed line"; then
    echo "✗ Test 5 failed: expected 'Warning: skipping invalid bed line' not found."
    echo "$OUT5"
    exit 1
fi

# Also, “chr1 50 40” => start > end => skip that region, but that shouldn't affect chr1:100
# The final VCF output should contain the line 'chr1 100'
if ! echo "$OUT5" | grep -q "^chr1[[:space:]]100"; then
    echo "✗ Test 5 failed: expected to see 'chr1\\t100' line in output"
    echo "$OUT5"
    exit 1
fi

echo "✓ Test 5 passed"

echo "✅ All tests for VCFX_region_subsampler passed successfully!"
