#!/usr/bin/env bash

# Exit immediately on any error
set -e
# Also fail if any subcommand in a pipe fails
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"

# Adjust if needed:
PARSER_BIN="${SCRIPT_DIR}/../build/src/VCFX_info_parser/VCFX_info_parser"

mkdir -p "${SCRIPT_DIR}/data/info_parser"

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
if ! "${PARSER_BIN}" --help 2>&1 | grep -q "VCFX_info_parser"; then
    echo "✗ Test failed: help - Did not find expected usage message."
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: Missing --info
###############################################################################
echo "Test 2: Missing --info"
if "${PARSER_BIN}" < /dev/null 2> "${SCRIPT_DIR}/data/info_parser/missing_info.err" || true; then
    :
fi

if ! grep -q "Error: INFO fields not specified" "${SCRIPT_DIR}/data/info_parser/missing_info.err"; then
    echo "✗ Test failed: missing --info - expected error not found"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Empty --info argument
###############################################################################
echo "Test 3: Empty --info"
if "${PARSER_BIN}" --info "" < /dev/null 2> "${SCRIPT_DIR}/data/info_parser/empty_info.err" || true; then
    :
fi

if ! grep -q "Error: INFO fields not specified" "${SCRIPT_DIR}/data/info_parser/empty_info.err"; then
    echo "✗ Test failed: empty --info - expected error not found"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Basic usage with a minimal VCF, single field
###############################################################################
echo "Test 4: Basic usage (single field)"
cat > "${SCRIPT_DIR}/data/info_parser/basic.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10
1	200	rs2	T	C	.	PASS	DP=20
EOF

"${PARSER_BIN}" --info "DP" < "${SCRIPT_DIR}/data/info_parser/basic.vcf" > "${SCRIPT_DIR}/data/info_parser/basic.out"

# Check header columns => CHROM POS ID REF ALT DP
if ! head -n1 "${SCRIPT_DIR}/data/info_parser/basic.out" | grep -q "^CHROM[[:space:]]POS[[:space:]]ID[[:space:]]REF[[:space:]]ALT[[:space:]]DP"; then
    echo "✗ Test failed: basic - header mismatch"
    exit 1
fi

# Check that we have the correct DP values
grep -q "^1[[:space:]]100[[:space:]].[[:space:]]A[[:space:]]G[[:space:]]10" "${SCRIPT_DIR}/data/info_parser/basic.out" \
  || (echo "✗ Test failed: basic - missing DP=10 line"; exit 1)

grep -q "^1[[:space:]]200[[:space:]]rs2[[:space:]]T[[:space:]]C[[:space:]]20" "${SCRIPT_DIR}/data/info_parser/basic.out" \
  || (echo "✗ Test failed: basic - missing DP=20 line"; exit 1)

echo "✓ Test 4 passed"

###############################################################################
# Test 5: Multiple fields (DP,AF)
###############################################################################
echo "Test 5: Multiple fields (DP,AF)"
cat > "${SCRIPT_DIR}/data/info_parser/multi.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10;AF=0.25
1	200	rs2	T	C	.	PASS	AF=0.75;DP=20
1	300	.	G	A	.	PASS	DP=15
EOF

"${PARSER_BIN}" --info "DP,AF" < "${SCRIPT_DIR}/data/info_parser/multi.vcf" > "${SCRIPT_DIR}/data/info_parser/multi.out"

# We expect columns: CHROM POS ID REF ALT DP AF
if ! head -n1 "${SCRIPT_DIR}/data/info_parser/multi.out" | grep -q "CHROM.*POS.*ID.*REF.*ALT.*DP.*AF"; then
    echo "✗ Test failed: multi - header mismatch"
    exit 1
fi

# lines
# 1) DP=10,AF=0.25
grep -q "^1[[:space:]]100[[:space:]].[[:space:]]A[[:space:]]G[[:space:]]10[[:space:]]0.25" "${SCRIPT_DIR}/data/info_parser/multi.out" \
  || (echo "✗ Test failed: multi - missing DP=10,AF=0.25"; exit 1)
# 2) AF=0.75,DP=20
grep -q "^1[[:space:]]200[[:space:]]rs2[[:space:]]T[[:space:]]C[[:space:]]20[[:space:]]0.75" "${SCRIPT_DIR}/data/info_parser/multi.out" \
  || (echo "✗ Test failed: multi - missing DP=20,AF=0.75"; exit 1)
# 3) DP=15 but no AF => AF => '.' 
grep -q "^1[[:space:]]300[[:space:]].[[:space:]]G[[:space:]]A[[:space:]]15[[:space:]]." "${SCRIPT_DIR}/data/info_parser/multi.out" \
  || (echo "✗ Test failed: multi - expected AF='.' for line 3"; exit 1)

echo "✓ Test 5 passed"

###############################################################################
# Test 6: Flags / missing fields
###############################################################################
echo "Test 6: Flags / missing fields"
cat > "${SCRIPT_DIR}/data/info_parser/flags.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10;SOMATIC;F2R1=xyz
1	200	rs2	C	T	.	PASS	AF=0.5
1	300	.	G	A	.	PASS	SOMATIC
EOF

"${PARSER_BIN}" --info "DP,AF,SOMATIC,F2R1" < "${SCRIPT_DIR}/data/info_parser/flags.vcf" > "${SCRIPT_DIR}/data/info_parser/flags.out"

# columns => CHROM POS ID REF ALT DP AF SOMATIC F2R1
head -n1 "${SCRIPT_DIR}/data/info_parser/flags.out" | grep -q "DP.*AF.*SOMATIC.*F2R1" \
  || (echo "✗ Test failed: flags - header mismatch"; exit 1)

# row1: DP=10, AF='.', SOMATIC has no '=', so it should appear as "", F2R1=xyz
grep -q "^1[[:space:]]100[[:space:]].[[:space:]]A[[:space:]]G[[:space:]]10[[:space:]].[[:space:]].[[:space:]]xyz" "${SCRIPT_DIR}/data/info_parser/flags.out" \
  || (echo "✗ Test failed: flags - row1 mismatch"; exit 1)

# row2: DP='.', AF=0.5, SOMATIC='.', F2R1='.'
grep -q "^1[[:space:]]200[[:space:]]rs2[[:space:]]C[[:space:]]T[[:space:]].[[:space:]]0.5[[:space:]].[[:space:]]." "${SCRIPT_DIR}/data/info_parser/flags.out" \
  || (echo "✗ Test failed: flags - row2 mismatch"; exit 1)

# row3: DP='.', AF='.', SOMATIC => '', F2R1='.'
grep -q "^1[[:space:]]300[[:space:]].[[:space:]]G[[:space:]]A[[:space:]].[[:space:]].[[:space:]].[[:space:]]." "${SCRIPT_DIR}/data/info_parser/flags.out" \
  || (echo "✗ Test failed: flags - row3 mismatch"; exit 1)

echo "✓ Test 6 passed"

###############################################################################
# Test 7: Invalid lines => parseInfoFields should warn + skip
###############################################################################
echo "Test 7: Invalid lines"
cat > "${SCRIPT_DIR}/data/info_parser/invalid.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10
ThisLineIsTotallyWrong
1	200	rs2	C	T	.	PASS	AF=0.5
EOF

"${PARSER_BIN}" --info "DP,AF" < "${SCRIPT_DIR}/data/info_parser/invalid.vcf" \
    2> "${SCRIPT_DIR}/data/info_parser/invalid.err" \
    > "${SCRIPT_DIR}/data/info_parser/invalid.out"

# We expect a warning for skipping invalid line
if ! grep -q "Warning: Skipping invalid VCF line" "${SCRIPT_DIR}/data/info_parser/invalid.err"; then
    echo "✗ Test failed: invalid lines - missing warning"
    exit 1
fi

# The 2 valid lines should appear in the output
grep -q "^1[[:space:]]100" "${SCRIPT_DIR}/data/info_parser/invalid.out" || (echo "✗ Test failed: invalid - missing line1"; exit 1)
grep -q "^1[[:space:]]200" "${SCRIPT_DIR}/data/info_parser/invalid.out" || (echo "✗ Test failed: invalid - missing line2"; exit 1)

echo "✓ Test 7 passed"

###############################################################################
# Test 8: CRLF line endings
###############################################################################
echo "Test 8: CRLF line endings"
cat > "${SCRIPT_DIR}/data/info_parser/crlf_unix.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=5
1	200	rs2	G	T	.	PASS	AF=0.25
EOF

sed 's/$/\r/' "${SCRIPT_DIR}/data/info_parser/crlf_unix.vcf" > "${SCRIPT_DIR}/data/info_parser/crlf.vcf"

"${PARSER_BIN}" --info "DP,AF" < "${SCRIPT_DIR}/data/info_parser/crlf.vcf" > "${SCRIPT_DIR}/data/info_parser/crlf.out"

# row1 => DP=5, AF='.'
grep -q "^1[[:space:]]100[[:space:]].[[:space:]]A[[:space:]]G[[:space:]]5[[:space:]]." "${SCRIPT_DIR}/data/info_parser/crlf.out" \
  || (echo "✗ Test failed: CRLF - missing row1 DP=5"; exit 1)
# row2 => DP='.', AF=0.25
grep -q "^1[[:space:]]200[[:space:]]rs2[[:space:]]G[[:space:]]T[[:space:]].[[:space:]]0.25" "${SCRIPT_DIR}/data/info_parser/crlf.out" \
  || (echo "✗ Test failed: CRLF - missing row2 AF=0.25"; exit 1)

echo "✓ Test 8 passed"

###############################################################################
# Test 9: Partial final line (no trailing newline)
###############################################################################
echo "Test 9: Partial final line"
cat > "${SCRIPT_DIR}/data/info_parser/partial_unix.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10
1	200	rs2	T	C	.	PASS	AF=0.3
EOF

truncate -s -1 "${SCRIPT_DIR}/data/info_parser/partial_unix.vcf"

"${PARSER_BIN}" --info "DP,AF" < "${SCRIPT_DIR}/data/info_parser/partial_unix.vcf" > "${SCRIPT_DIR}/data/info_parser/partial_unix.out"

# row1 => DP=10, AF='.'
grep -q "^1[[:space:]]100[[:space:]].[[:space:]]A[[:space:]]G[[:space:]]10[[:space:]]." "${SCRIPT_DIR}/data/info_parser/partial_unix.out" \
  || (echo "✗ Test failed: partial line - row1 mismatch"; exit 1)

# row2 => DP='.', AF=0.3
grep -q "^1[[:space:]]200[[:space:]]rs2[[:space:]]T[[:space:]]C[[:space:]].[[:space:]]0.3" "${SCRIPT_DIR}/data/info_parser/partial_unix.out" \
  || (echo "✗ Test failed: partial line - row2 mismatch"; exit 1)

echo "✓ Test 9 passed"

###############################################################################
echo "✓ All tests passed"
