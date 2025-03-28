#!/usr/bin/env bash

set -e
set -o pipefail
# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Location of the built indexer binary (adjust if your build directory differs)
CALCULATOR_BIN="${SCRIPT_DIR}/../build/src/VCFX_ld_calculator/VCFX_ld_calculator"

###############################################################################
# Test 1: Help/usage
###############################################################################
echo "Test 1: Help / usage"
if ! "$CALCULATOR_BIN" --help 2>&1 | grep -q "VCFX_ld_calculator: Calculate pairwise LD"; then
    echo "✗ Test failed: help - usage message not found."
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: Single variant => no pairwise LD
###############################################################################
echo "Test 2: Single variant"

cat > single.vcf <<EOF
##fileformat=VCFv4.2
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
1       100     .       A       G       .       PASS    .       GT      0/0
EOF

# No --region => use entire file
OUT2=$(cat single.vcf | "$CALCULATOR_BIN")
# We expect #LD_MATRIX_START plus message about "No or only one variant"
if ! echo "$OUT2" | grep -q "No or only one variant in the region => no pairwise LD."; then
    echo "✗ Test failed: single variant => expected 'no pairwise LD' message not found."
    echo "Output was:"
    echo "$OUT2"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Two variants => produce 2x2 matrix
###############################################################################
echo "Test 3: Two variants"

cat > two_variants.vcf <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t200\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\t1/1
EOF

OUT3=$(cat two_variants.vcf | "$CALCULATOR_BIN")

# We expect:
# - two data lines
# - #LD_MATRIX_START
# - 2x2 matrix with "1:100" and "1:200"
if ! echo "$OUT3" | grep -q "#LD_MATRIX_START"; then
    echo "✗ Test failed: two variants => missing #LD_MATRIX_START"
    echo "Output was:"
    echo "$OUT3"
    exit 1
fi
if ! echo "$OUT3" | grep -q "^1:100\s"; then
    echo "✗ Test failed: two variants => missing row for 1:100"
    echo "Output was:"
    echo "$OUT3"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Region excludes everything => no pairwise LD
###############################################################################
echo "Test 4: Region excludes everything"

cat > out_of_range.vcf <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0
1\t200\t.\tT\tC\t.\tPASS\t.\tGT\t0/1
EOF

# region is chr1:300-400 => excludes POS=100 & 200 => no variants in region
OUT4=$(cat out_of_range.vcf | "$CALCULATOR_BIN" --region 1:300-400)

if ! echo "$OUT4" | grep -q "No or only one variant in the region => no pairwise LD."; then
    echo "✗ Test failed: region excludes everything => expected 'no pairwise LD' message not found."
    echo "Output was:"
    echo "$OUT4"
    exit 1
fi
echo "✓ Test 4 passed"

###############################################################################
# Test 5: Region includes a subset
###############################################################################
echo "Test 5: Region includes a subset"

cat > partial_region.vcf <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
1\t250\t.\tT\tC\t.\tPASS\t.\tGT\t1/1
1\t350\t.\tG\tT\t.\tPASS\t.\tGT\t0/0
EOF

# We'll only keep variants at [200..300], so that includes POS=250 => 1 variant
OUT5=$(cat partial_region.vcf | "$CALCULATOR_BIN" --region 1:200-300)

if ! echo "$OUT5" | grep -q "No or only one variant in the region => no pairwise LD."; then
    echo "✗ Test failed: region subset => expected 'no pairwise LD' message"
    echo "Output was:"
    echo "$OUT5"
    exit 1
fi
echo "✓ Test 5 passed"

###############################################################################
# Test 6: Missing / multi-allelic genotypes
###############################################################################
echo "Test 6: Missing or multi-allelic"

cat > missing_multi.vcf <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t100\tvarA\tA\tG\t.\tPASS\t.\tGT\t0/1\t./.\t1/1
1\t200\tvarB\tT\tC\t.\tPASS\t.\tGT\t1/2\t0/1\t0/0
1\t300\tvarC\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0
EOF

OUT6=$(cat missing_multi.vcf | "$CALCULATOR_BIN")

# We expect #LD_MATRIX_START and a 3x3 matrix
if ! echo "$OUT6" | grep -q "#LD_MATRIX_START"; then
    echo "✗ Test failed: missing/multi => no #LD_MATRIX_START found"
    echo "$OUT6"
    exit 1
fi
if ! echo "$OUT6" | grep -q "1:200"; then
    echo "✗ Test failed: missing/multi => variant 1:200 not found in LD matrix"
    echo "$OUT6"
    exit 1
fi
echo "✓ Test 6 passed"

###############################################################################
# Test 7: Data line before #CHROM => expect error message
###############################################################################
echo "Test 7: Data line before #CHROM"

cat > no_header.vcf <<EOF
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t200\t.\tA\tC\t.\tPASS\t.\tGT\t0/1
EOF

ERR7=$(cat no_header.vcf | "$CALCULATOR_BIN" 2>&1 || true)

# The code prints "Error: encountered data line before #CHROM."
if ! echo "$ERR7" | grep -q "Error: encountered data line before #CHROM."; then
    echo "✗ Test failed: data line before #CHROM => expected error message not found."
    echo "Stderr was:"
    echo "$ERR7"
    exit 1
fi
echo "✓ Test 7 passed"
