#!/usr/bin/env bash

# Exit immediately on error
set -e
# Fail if any subcommand in a pipe fails
set -o pipefail

# Directory of this script
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"

# Location of the built aggregator binary (adjust if needed)
AGG_BIN="${SCRIPT_DIR}/../build/src/VCFX_info_aggregator/VCFX_info_aggregator"

# Ensure we have a subdirectory for aggregator tests
mkdir -p "${SCRIPT_DIR}/data/aggregator"

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
if ! "${AGG_BIN}" --help 2>&1 | grep -q "VCFX_info_aggregator"; then
    echo "✗ Test failed: help - Did not find expected usage message."
    exit 1
fi
echo "✓ Test 1 passed"


###############################################################################
# Test 2: Missing --aggregate-info
###############################################################################
echo "Test 2: Missing --aggregate-info"
# We expect an error message
if "${AGG_BIN}" < /dev/null 2> "${SCRIPT_DIR}/data/aggregator/missing_arg.err" || true; then
    : # aggregator might return non-zero
fi
if ! grep -q "Error: Must specify --aggregate-info" "${SCRIPT_DIR}/data/aggregator/missing_arg.err"; then
    echo "✗ Test failed: missing --aggregate-info - did not find expected error"
    exit 1
fi
echo "✓ Test 2 passed"


###############################################################################
# Test 3: Empty --aggregate-info
###############################################################################
echo "Test 3: Empty --aggregate-info"
# We give an empty argument, aggregator should fail with "Error: no valid fields"
if "${AGG_BIN}" --aggregate-info "" < /dev/null 2> "${SCRIPT_DIR}/data/aggregator/empty_fields.err" || true; then
    :
fi
if ! grep -q "Error: Must specify --aggregate-info with at least one field." \
       "${SCRIPT_DIR}/data/aggregator/empty_fields.err"; then
    echo "✗ Test failed: empty --aggregate-info - missing error"
    exit 1
fi
echo "✓ Test 3 passed"


###############################################################################
# Test 4: Basic aggregator with a single numeric field (DP)
###############################################################################
echo "Test 4: Basic aggregator (DP)"
cat > "${SCRIPT_DIR}/data/aggregator/basic.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10
1	200	.	T	C	.	PASS	DP=20
EOF

"${AGG_BIN}" --aggregate-info "DP" < "${SCRIPT_DIR}/data/aggregator/basic.vcf" > "${SCRIPT_DIR}/data/aggregator/basic.out"

# Verify original lines are present
if ! grep -q "^1	100" "${SCRIPT_DIR}/data/aggregator/basic.out"; then
    echo "✗ Test failed: basic aggregator - missing original line"
    exit 1
fi

# Verify aggregator summary
if ! grep -q "^#AGGREGATION_SUMMARY" "${SCRIPT_DIR}/data/aggregator/basic.out"; then
    echo "✗ Test failed: basic aggregator - missing summary"
    exit 1
fi

# Check for DP summary. We expect sum=30, average=15
if ! grep -q "^DP: Sum=30, Average=15" "${SCRIPT_DIR}/data/aggregator/basic.out"; then
    echo "✗ Test failed: basic aggregator - wrong DP sum/average"
    exit 1
fi
echo "✓ Test 4 passed"


###############################################################################
# Test 5: Multiple numeric fields (DP,AF)
###############################################################################
echo "Test 5: Multiple numeric fields (DP,AF)"
cat > "${SCRIPT_DIR}/data/aggregator/multi.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10;AF=0.25
1	200	.	T	C	.	PASS	DP=20;AF=0.75
1	300	.	G	A	.	PASS	DP=5;AF=0.10
EOF

"${AGG_BIN}" --aggregate-info "DP,AF" < "${SCRIPT_DIR}/data/aggregator/multi.vcf" > "${SCRIPT_DIR}/data/aggregator/multi.out"

# Summaries:
#   DP: sum=35 (10+20+5), average=11.6667
#   AF: sum=1.10 (0.25+0.75+0.10), average=0.3667
# We'll just check sums, and approximate average to 4 decimal places.

if ! grep -q "#AGGREGATION_SUMMARY" "${SCRIPT_DIR}/data/aggregator/multi.out"; then
    echo "✗ Test failed: multi aggregator - no summary"
    exit 1
fi

# Check DP
if ! grep -q "^DP: Sum=35," "${SCRIPT_DIR}/data/aggregator/multi.out"; then
    echo "✗ Test failed: multi aggregator - DP sum not 35"
    exit 1
fi

# Check AF
if ! grep -q "^AF: Sum=1.1," "${SCRIPT_DIR}/data/aggregator/multi.out"; then
    echo "✗ Test failed: multi aggregator - AF sum not 1.1"
    exit 1
fi
echo "✓ Test 5 passed"


###############################################################################
# Test 6: Non-numeric values (skip them)
###############################################################################
echo "Test 6: Non-numeric or partial numeric"
cat > "${SCRIPT_DIR}/data/aggregator/non_numeric.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=abc;AF=0.2
1	200	.	C	T	.	PASS	DP=20;AF=NaN
1	300	.	G	A	.	PASS	DP=30;AF=0.3
EOF

"${AGG_BIN}" --aggregate-info "DP,AF" < "${SCRIPT_DIR}/data/aggregator/non_numeric.vcf" > "${SCRIPT_DIR}/data/aggregator/non_numeric.out"

# We skip "DP=abc" and "AF=NaN" as they're not parseable
# DP: sum=50 (20+30), avg=25
# AF: sum=0.5 (0.2+0.3), avg=0.25

if ! grep -q "DP: Sum=50," "${SCRIPT_DIR}/data/aggregator/non_numeric.out"; then
    echo "✗ Test failed: non_numeric - expected DP sum=50"
    exit 1
fi
if ! grep -q "AF: Sum=0.5," "${SCRIPT_DIR}/data/aggregator/non_numeric.out"; then
    echo "✗ Test failed: non_numeric - expected AF sum=0.5"
    exit 1
fi
echo "✓ Test 6 passed"


###############################################################################
# Test 7: No #CHROM header => aggregator should error
###############################################################################
echo "Test 7: No #CHROM header"
cat > "${SCRIPT_DIR}/data/aggregator/no_chrom.vcf" << EOF
1	100	rs1	A	G	.	PASS	DP=10
1	200	rs2	C	T	.	PASS	DP=20
EOF

# aggregator should detect data lines before #CHROM and print error
if "${AGG_BIN}" --aggregate-info "DP" < "${SCRIPT_DIR}/data/aggregator/no_chrom.vcf" 2> "${SCRIPT_DIR}/data/aggregator/no_chrom.err" || true; then
    :
fi

if ! grep -q "Error: encountered data line before #CHROM header" "${SCRIPT_DIR}/data/aggregator/no_chrom.err"; then
    echo "✗ Test failed: no_chrom - missing expected error"
    exit 1
fi
echo "✓ Test 7 passed"


###############################################################################
# Test 8: CRLF line endings
###############################################################################
echo "Test 8: CRLF line endings"
cat > "${SCRIPT_DIR}/data/aggregator/crlf_unix.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=5;AF=0.10
1	200	.	T	C	.	PASS	DP=15;AF=0.30
EOF

# Convert to CRLF by appending \r
sed 's/$/\r/' "${SCRIPT_DIR}/data/aggregator/crlf_unix.vcf" > "${SCRIPT_DIR}/data/aggregator/crlf.vcf"

"${AGG_BIN}" --aggregate-info "DP,AF" < "${SCRIPT_DIR}/data/aggregator/crlf.vcf" > "${SCRIPT_DIR}/data/aggregator/crlf.out"

# sums: DP=20, AF=0.4
if ! grep -q "#AGGREGATION_SUMMARY" "${SCRIPT_DIR}/data/aggregator/crlf.out"; then
    echo "✗ Test failed: crlf - no summary"
    exit 1
fi
if ! grep -q "DP: Sum=20," "${SCRIPT_DIR}/data/aggregator/crlf.out"; then
    echo "✗ Test failed: crlf - sum=20 not found"
    exit 1
fi
if ! grep -q "AF: Sum=0.4," "${SCRIPT_DIR}/data/aggregator/crlf.out"; then
    echo "✗ Test failed: crlf - AF sum=0.4 not found"
    exit 1
fi
echo "✓ Test 8 passed"


###############################################################################
# Test 9: Partial final line (no trailing newline)
###############################################################################
echo "Test 9: Partial final line with no trailing newline"
cat > "${SCRIPT_DIR}/data/aggregator/partial_unix.vcf" << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10;AF=0.2
1	200	.	C	T	.	PASS	DP=40;AF=0.3
EOF

# remove final newline
truncate -s -1 "${SCRIPT_DIR}/data/aggregator/partial_unix.vcf"

"${AGG_BIN}" --aggregate-info "DP,AF" < "${SCRIPT_DIR}/data/aggregator/partial_unix.vcf" > "${SCRIPT_DIR}/data/aggregator/partial_unix.out"

# sums: DP=50, AF=0.5
if ! grep -q "DP: Sum=50," "${SCRIPT_DIR}/data/aggregator/partial_unix.out"; then
    echo "✗ Test failed: partial line - DP sum=50 not found"
    exit 1
fi
if ! grep -q "AF: Sum=0.5," "${SCRIPT_DIR}/data/aggregator/partial_unix.out"; then
    echo "✗ Test failed: partial line - AF sum=0.5 not found"
    exit 1
fi
echo "✓ Test 9 passed"

###############################################################################
echo "✓ All tests passed"
