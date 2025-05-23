#!/usr/bin/env bash

# Exit immediately on any error
set -e
# Also fail if any command in a pipe fails
set -o pipefail
# Print each command in debug mode (optional, comment out if too verbose):
#set -x

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
INDEXER_BIN="${SCRIPT_DIR}/../build/src/VCFX_indexer/VCFX_indexer"

# Ensure data directory
mkdir -p "${SCRIPT_DIR}/data/indexer"

# A helper function that writes only space‐separated lines (no tabs).
write_space_only() {
    local dest="$1"
    shift
    # Replace any tabs with spaces, write to $dest
    sed 's/\t/ /g' > "$dest"

    # Verify no tabs remain
    if grep -q $'\t' "$dest"; then
        echo "✗ Test failed: found unexpected tab in $dest!"
        exit 1
    fi
}

###############################################################################
# Test 1: Basic indexing
###############################################################################
echo "Test 1: Basic indexing"
cat > "${SCRIPT_DIR}/data/indexer/basic.vcf" <<EOF
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	100	PASS	AC=2;AN=4	GT	0/0	1/1
1	200	rs2	C	T	100	PASS	AC=2;AN=4	GT	0/1	0/1
1	300	rs3	G	A	100	PASS	AC=2;AN=4	GT	0/0	1/1
EOF

# Run the indexer
"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/basic.vcf" > "${SCRIPT_DIR}/data/indexer/basic.index"

# Check the header
if ! grep -q "^CHROM[[:space:]]POS[[:space:]]FILE_OFFSET" "${SCRIPT_DIR}/data/indexer/basic.index"; then
    echo "✗ Test failed: basic - Missing header line"
    exit 1
fi

# Expect 3 lines of data
INDEX_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/basic.index")
if [ "$INDEX_COUNT" -ne 3 ]; then
    echo "✗ Test failed: basic - Expected 3 entries, got $INDEX_COUNT"
    exit 1
fi

echo "✓ Test 1 passed"

###############################################################################
# Test 2: FILE_OFFSET accuracy
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
# Test 3: Missing #CHROM header
###############################################################################
echo "Test 3: Missing #CHROM header"
cat > "${SCRIPT_DIR}/data/indexer/no_header.vcf" <<EOF
1	100	rs1	A	G	100	PASS	AC=2;AN=4	GT	0/0	1/1
1	200	rs2	C	T	100	PASS	AC=2;AN=4	GT	0/1	0/1
EOF

"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/no_header.vcf" 2> "${SCRIPT_DIR}/data/indexer/no_header.err" || true

if ! grep -q "Error: no #CHROM header found" "${SCRIPT_DIR}/data/indexer/no_header.err"; then
    echo "✗ Test failed: no_header - Expected error message not found"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Empty lines and comments
###############################################################################
echo "Test 4: Empty lines and comments"
cat > "${SCRIPT_DIR}/data/indexer/with_comments.vcf" <<EOF
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2

1	100	rs1	A	G	100	PASS	AC=2;AN=4	GT	0/0	1/1

# This is a comment
1	200	rs2	C	T	100	PASS	AC=2;AN=4	GT	0/1	0/1

EOF

"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/with_comments.vcf" > "${SCRIPT_DIR}/data/indexer/with_comments.index"
INDEX_COUNT=$(grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/with_comments.index")
if [ "$INDEX_COUNT" -ne 2 ]; then
    echo "✗ Test failed: with_comments - Expected 2 entries, got $INDEX_COUNT"
    exit 1
fi
echo "✓ Test 4 passed"

###############################################################################
# Test 5: Non-seekable input (pipe)
###############################################################################
echo "Test 5: Non-seekable input"
cat "${SCRIPT_DIR}/data/indexer/basic.vcf" | "${INDEXER_BIN}" > "${SCRIPT_DIR}/data/indexer/pipe.index"

if ! diff "${SCRIPT_DIR}/data/indexer/basic.index" "${SCRIPT_DIR}/data/indexer/pipe.index" > /dev/null; then
    echo "✗ Test failed: pipe - Index differs from file input"
    exit 1
fi
echo "✓ Test 5 passed"

###############################################################################
# Test 6: Malformed VCF (space-separated, no tabs)
###############################################################################
echo "Test 6: Malformed VCF"
write_space_only "${SCRIPT_DIR}/data/indexer/malformed.vcf" <<EOF
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROM  POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   SAMPLE1   SAMPLE2
1  100  rs1   A    G     100   PASS    AC=2;AN=4   GT   0/0   1/1
1  200  rs2   C    T     100   PASS    AC=2;AN=4   GT   0/1   0/1
EOF

"${INDEXER_BIN}" < "${SCRIPT_DIR}/data/indexer/malformed.vcf" > "${SCRIPT_DIR}/data/indexer/malformed.index"

INDEX_COUNT=$(
  grep -c "^1[[:space:]]" "${SCRIPT_DIR}/data/indexer/malformed.index" || true
)

if [ "$INDEX_COUNT" -ne 0 ]; then
    echo "✗ Test failed: malformed - Expected 0 index entries, got $INDEX_COUNT"
    exit 1
fi
echo "✓ Test 6 passed"

###############################################################################
# Test 7: CRLF line endings
###############################################################################
echo "Test 7: CRLF line endings"
cat > "${SCRIPT_DIR}/data/indexer/crlf_unix.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	rs1	A	G	100	PASS	.	GT	0/1
1	200	rs2	T	C	100	PASS	.	GT	0/0
EOF

sed 's/$/\r/' "${SCRIPT_DIR}/data/indexer/crlf_unix.vcf" > "${SCRIPT_DIR}/data/indexer/crlf.vcf"

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
# Test 8: Partial final line
###############################################################################
echo "Test 8: Partial final line with no trailing newline"
cat > "${SCRIPT_DIR}/data/indexer/partial_unix.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	rs1	A	G	100	PASS	.	GT	0/1
1	200	rs2	T	C	100	PASS	.	GT	0/0
EOF

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

echo "✓ All tests passed"
