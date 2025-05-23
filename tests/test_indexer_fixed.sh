#!/usr/bin/env bash

# Exit immediately on any error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Location of the built indexer binary (adjust if your build directory differs)
INDEXER_BIN="${SCRIPT_DIR}/../build/src/VCFX_indexer/VCFX_indexer"

# Ensure data directory exists
mkdir -p "${SCRIPT_DIR}/data/indexer"

function write_space_only() {
    local dest="$1"
    # Shift off the filename, so any 'here-doc' lines come from stdin
    shift

    # 1) Replace tabs with spaces
    # 2) Write to $dest
    sed 's/\t/ /g' > "$dest"

    # 3) Double check that no tabs remain
    if grep -q $'\t' "$dest"; then
        echo "✗ Test failed: found unexpected tab in $dest!"
        exit 1
    fi
}

###############################################################################
# Test 1: Basic indexing with a properly formatted VCF file
###############################################################################
echo "Test 1: Basic indexing"
cat > "${SCRIPT_DIR}/data/indexer/basic.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROMPOSIDREFALTQUALFILTERINFOFORMATSAMPLE1SAMPLE2
1100rs1AG100PASSAC=2;AN=4GT0/01/1
1200rs2CT100PASSAC=2;AN=4GT0/10/1
1300rs3GA100PASSAC=2;AN=4GT0/01/1
EOF

# Run the indexer
"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/basic.vcf" > "${SCRIPT_DIR}/data/indexer/basic.index"

# Check the header
if ! grep -q "^CHROM[[:space:]]POS[[:space:]]FILE_OFFSET" "${SCRIPT_DIR}/data/indexer/basic.index"; then
    echo "✗ Test failed: basic - Missing header line"
    exit 1
fi

# Count variant lines
INDEX_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/basic.index")
if [ "$INDEX_COUNT" -ne 3 ]; then
    echo "✗ Test failed: basic - Expected 3 index entries, got $INDEX_COUNT"
    exit 1
fi

echo "✓ Test 1 passed"


###############################################################################
# Test 2: Verify FILE_OFFSET is accurate (using dd to seek raw file bytes)
###############################################################################
echo "Test 2: FILE_OFFSET accuracy"
FIRST_OFFSET=$(grep "^1[[:space:]]100" "${SCRIPT_DIR}/data/indexer/basic.index" | cut -f3)

FIRST_LINE=$(
  dd if="${SCRIPT_DIR}/data/indexer/basic.vcf" bs=1 skip="$FIRST_OFFSET" count=100 2>/dev/null \
    | tr -d '\r\n' \
    | sed 's/\t/ /g' \
    | sed 's/  */ /g' \
    | sed 's/^ //;s/ $//' \
    | cut -d' ' -f1-3
)

if [ "$FIRST_LINE" != "1 100 rs1" ]; then
    echo "✗ Test failed: basic - FILE_OFFSET is incorrect"
    echo "Expected '1 100 rs1', got: '$FIRST_LINE'"
    exit 1
fi

echo "✓ Test 2 passed"


###############################################################################
# Test 3: Handle missing #CHROM header
###############################################################################
echo "Test 3: Missing #CHROM header"
cat > "${SCRIPT_DIR}/data/indexer/no_header.vcf" << 'EOF'
1100rs1AG100PASSAC=2;AN=4GT0/01/1
1200rs2CT100PASSAC=2;AN=4GT0/10/1
EOF

# Capture stderr to detect error message
"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/no_header.vcf" 2> "${SCRIPT_DIR}/data/indexer/no_header.err" || true

if ! grep -q "Error: no #CHROM header found" "${SCRIPT_DIR}/data/indexer/no_header.err"; then
    echo "✗ Test failed: no_header - Expected error message not found"
    exit 1
fi

echo "✓ Test 3 passed"


###############################################################################
# Test 4: Handle empty lines and comments
###############################################################################
echo "Test 4: Empty lines and comments"
cat > "${SCRIPT_DIR}/data/indexer/with_comments.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROMPOSIDREFALTQUALFILTERINFOFORMATSAMPLE1SAMPLE2

1100rs1AG100PASSAC=2;AN=4GT0/01/1

# This is a comment
1200rs2CT100PASSAC=2;AN=4GT0/10/1

EOF

"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/with_comments.vcf" > "${SCRIPT_DIR}/data/indexer/with_comments.index"

INDEX_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/with_comments.index")
if [ "$INDEX_COUNT" -ne 2 ]; then
    echo "✗ Test failed: with_comments - Expected 2 index entries, got $INDEX_COUNT"
    exit 1
fi

echo "✓ Test 4 passed"


###############################################################################
# Test 5: Handle non-seekable input (pipe)
###############################################################################
echo "Test 5: Non-seekable input"
cat "${SCRIPT_DIR}/data/indexer/basic.vcf" | "${INDEXER_BIN}" > "${SCRIPT_DIR}/data/indexer/pipe.index"

if ! diff "${SCRIPT_DIR}/data/indexer/basic.index" "${SCRIPT_DIR}/data/indexer/pipe.index" > /dev/null; then
    echo "✗ Test failed: pipe - Index differs from file input"
    exit 1
fi

echo "✓ Test 5 passed"

###############################################################################
# Test 6: Handle malformed VCF (space-separated instead of tab-separated)
###############################################################################
echo "Test 6: Malformed VCF"

# Use our write_space_only function to forcibly remove any accidental tabs:
write_space_only "${SCRIPT_DIR}/data/indexer/malformed.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROM  POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   SAMPLE1   SAMPLE2
1  100  rs1   A    G     100   PASS    AC=2;AN=4   GT   0/0   1/1
1  200  rs2   C    T     100   PASS    AC=2;AN=4   GT   0/1   0/1
EOF

# Now run the indexer
"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/malformed.vcf" > "${SCRIPT_DIR}/data/indexer/malformed.index"

# We expect 0 lines starting with "1", meaning 0 variants
INDEX_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/malformed.index")
if [ "$INDEX_COUNT" -ne 0 ]; then
    echo "✗ Test failed: malformed - Expected 0 index entries, got $INDEX_COUNT"
    exit 1
fi

echo "✓ Test 6 passed"

###############################################################################
# Test 7: CRLF line endings
# We will transform a normal VCF into CRLF and check offsets
###############################################################################
echo "Test 7: CRLF line endings"

cat > "${SCRIPT_DIR}/data/indexer/crlf_unix.vcf" << 'EOF'
##fileformat=VCFv4.2
#CHROMPOSIDREFALTQUALFILTERINFOFORMATSAMPLE
1100rs1AG100PASS.GT0/1
1200rs2TC100PASS.GT0/0
EOF

# Convert to CRLF by appending \r before \n
sed 's/$/\r/' "${SCRIPT_DIR}/data/indexer/crlf_unix.vcf" > "${SCRIPT_DIR}/data/indexer/crlf.vcf"

# Index the CRLF file
"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/crlf.vcf" > "${SCRIPT_DIR}/data/indexer/crlf.index"

CRLF_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/crlf.index")
if [ "$CRLF_COUNT" -ne 2 ]; then
    echo "✗ Test failed: CRLF - Expected 2 index entries, got $CRLF_COUNT"
    exit 1
fi

CRLF_FIRST_OFFSET=$(grep "^1[[:space:]]100" "${SCRIPT_DIR}/data/indexer/crlf.index" | cut -f3)
CRLF_FIRST_LINE=$(
  dd if="${SCRIPT_DIR}/data/indexer/crlf.vcf" bs=1 skip="$CRLF_FIRST_OFFSET" count=50 2>/dev/null \
    | tr -d '\r\n' \
    | sed 's/\t/ /g' \
    | sed 's/  */ /g' \
    | sed 's/^ //;s/ $//' \
    | cut -d' ' -f1-3
)
if [ "$CRLF_FIRST_LINE" != "1 100 rs1" ]; then
    echo "✗ Test failed: CRLF - FILE_OFFSET is incorrect for the first variant"
    echo "Expected '1\\t100\\trs1', got: '$CRLF_FIRST_LINE'"
    exit 1
fi

echo "✓ Test 7 passed"


###############################################################################
# Test 8: Partial final line (no trailing newline)
###############################################################################
echo "Test 8: Partial final line with no trailing newline"

cat > "${SCRIPT_DIR}/data/indexer/partial_unix.vcf" << 'EOF'
##fileformat=VCFv4.2
#CHROMPOSIDREFALTQUALFILTERINFOFORMATSAMPLE
1100rs1AG100PASS.GT0/1
1200rs2TC100PASS.GT0/0
EOF

# Remove the final newline so the last variant line is partial
truncate -s -1 "${SCRIPT_DIR}/data/indexer/partial_unix.vcf"

"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/partial_unix.vcf" > "${SCRIPT_DIR}/data/indexer/partial_unix.index"

PARTIAL_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/partial_unix.index")
if [ "$PARTIAL_COUNT" -ne 2 ]; then
    echo "✗ Test failed: partial line - Expected 2 index entries, got $PARTIAL_COUNT"
    exit 1
fi

PARTIAL_OFFSET=$(grep "^1[[:space:]]200" "${SCRIPT_DIR}/data/indexer/partial_unix.index" | cut -f3)
PARTIAL_LINE=$(
  dd if="${SCRIPT_DIR}/data/indexer/partial_unix.vcf" bs=1 skip="$PARTIAL_OFFSET" count=50 2>/dev/null \
    | tr -d '\r\n' \
    | sed 's/\t/ /g' \
    | sed 's/  */ /g' \
    | sed 's/^ //;s/ $//' \
    | cut -d' ' -f1-3
)
if [ "$PARTIAL_LINE" != "1 200 rs2" ]; then
    echo "✗ Test failed: partial line - FILE_OFFSET is incorrect"
    echo "Expected '1\\t200\\trs2', got: '$PARTIAL_LINE'"
    exit 1
fi

echo "✓ Test 8 passed"

###############################################################################
echo "✓ All tests passed"
