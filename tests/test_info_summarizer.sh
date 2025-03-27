#!/usr/bin/env bash

# Exit immediately on any error
set -e
# Fail if a subcommand in a pipe fails
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
SUMMARIZER_BIN="${SCRIPT_DIR}/../build/src/VCFX_info_summarizer/VCFX_info_summarizer"

mkdir -p "${SCRIPT_DIR}/data/info_summarizer"

###############################################################################
# Utility: debugHex function
# This prints a debug hexdump of the given file so we can see if tabs/spaces
# are correct
###############################################################################
debugHex() {
    local file="$1"
    echo "[DEBUG] od -c of '${file}':"
    od -c "${file}" || true
    echo ""
}

###############################################################################
# Test 1: Help / usage
###############################################################################
echo "Test 1: Help / usage"
if ! "${SUMMARIZER_BIN}" --help 2>&1 | grep -q "VCFX_info_summarizer"; then
    echo "✗ Test failed: help - usage message not found."
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: Missing --info
###############################################################################
echo "Test 2: Missing --info"
if "${SUMMARIZER_BIN}" < /dev/null 2> "${SCRIPT_DIR}/data/info_summarizer/missing_info.err" || true; then
    :
fi

if ! grep -q "Error: INFO fields not specified" "${SCRIPT_DIR}/data/info_summarizer/missing_info.err"; then
    echo "✗ Test failed: missing --info - expected error not found"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Empty --info argument
###############################################################################
echo "Test 3: Empty --info"
if "${SUMMARIZER_BIN}" --info "" < /dev/null 2> "${SCRIPT_DIR}/data/info_summarizer/empty_info.err" || true; then
    :
fi

if ! grep -q "Error: INFO fields not specified" "${SCRIPT_DIR}/data/info_summarizer/empty_info.err"; then
    echo "✗ Test failed: empty --info - expected error not found"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Basic scenario with a single numeric field (DP)
###############################################################################
echo "Test 4: Basic scenario (DP)"

# We’ll create "basic.vcf" with actual tab characters (\t).
{
  echo "##fileformat=VCFv4.2"
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  echo -e "1\t100\t.\tA\tG\t.\tPASS\tDP=10"
  echo -e "1\t200\trs2\tT\tC\t.\tPASS\tDP=20"
  echo -e "1\t300\t.\tG\tA\t.\tPASS\tDP=30"
} > "${SCRIPT_DIR}/data/info_summarizer/basic.vcf"

# Debug: show a hexdump so we can confirm the file has tabs.
debugHex "${SCRIPT_DIR}/data/info_summarizer/basic.vcf"

# Now run the summarizer with "--info DP"
"${SUMMARIZER_BIN}" --info "DP" \
  < "${SCRIPT_DIR}/data/info_summarizer/basic.vcf" \
  > "${SCRIPT_DIR}/data/info_summarizer/basic.out"

# We expect DP => [10,20,30] => mean=20, median=20, mode=10
# Relax the grep pattern to allow either "20" or "20.0000", etc.
if ! grep -Eq "^DP[[:space:]]+(20(\.0000)?)[[:space:]]+(20(\.0000)?)[[:space:]]+(10(\.0000)?)" \
    "${SCRIPT_DIR}/data/info_summarizer/basic.out"; then
    
    echo "✗ Test failed: basic - DP stats mismatch"
    echo "Expected mean=20, median=20, mode=10 (with optional .0000)."
    echo "[DEBUG] Content of basic.out:"
    cat "${SCRIPT_DIR}/data/info_summarizer/basic.out"
    exit 1
fi

echo "✓ Test 4 passed"


###############################################################################
# Test 5: Multiple fields (DP,AF)
###############################################################################
echo "Test 5: Multiple fields (DP,AF)"
cat > "${SCRIPT_DIR}/data/info_summarizer/multi.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10;AF=0.25
1	200	rs2	C	T	.	PASS	DP=10;AF=0.75
1	300	.	G	A	.	PASS	DP=20;AF=0.75
1	400	.	A	G	.	PASS	DP=40
EOF

debugHex "${SCRIPT_DIR}/data/info_summarizer/multi.vcf"

"${SUMMARIZER_BIN}" --info "DP,AF" < "${SCRIPT_DIR}/data/info_summarizer/multi.vcf" > "${SCRIPT_DIR}/data/info_summarizer/multi.out"

# DP => [10, 10, 20, 40] => mean=20, median=15, mode=10
# AF => [0.25, 0.75, 0.75] => mean~0.5833, median=0.75, mode=0.75

# Relaxed pattern for DP: allow either `20` or `20.0000`, etc.
if ! grep -Eq "^DP[[:space:]]+(20(\.0000)?)[[:space:]]+(15(\.0000)?)[[:space:]]+(10(\.0000)?)" \
    "${SCRIPT_DIR}/data/info_summarizer/multi.out"; then
    echo "✗ Test failed: multi - DP stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/multi.out"
    exit 1
fi

# Relaxed pattern for AF:
# We want mean=0.5833, median=0.75, mode=0.75, but allow `.0000` suffix on the 0.75 fields.
if ! grep -Eq "^AF[[:space:]]*0\.5833[[:space:]]+(0\.75(\.0000)?)[[:space:]]+(0\.75(\.0000)?)" \
    "${SCRIPT_DIR}/data/info_summarizer/multi.out"; then
    echo "✗ Test failed: multi - AF stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/multi.out"
    exit 1
fi

echo "✓ Test 5 passed"

###############################################################################
# Test 6: Non-numeric, multiple comma-separated values
###############################################################################
echo "Test 6: Non-numeric and multi-values"
cat > "${SCRIPT_DIR}/data/info_summarizer/non_numeric.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10,abc;AF=0.5
1	200	.	C	T	.	PASS	DP=20,30;AF=NaN
1	300	.	G	A	.	PASS	DP=40;AF=0.2,0.3
EOF

debugHex "${SCRIPT_DIR}/data/info_summarizer/non_numeric.vcf"

"${SUMMARIZER_BIN}" --info "DP,AF" < "${SCRIPT_DIR}/data/info_summarizer/non_numeric.vcf" > "${SCRIPT_DIR}/data/info_summarizer/non_numeric.out"

# DP => [10,20,30,40] => mean=25, median=25, mode=10
# AF => [0.5,0.2,0.3] => mean=0.3333, median=0.3, mode=0.5
if ! grep -q "^DP[[:space:]]*25.0000[[:space:]]*25[[:space:]]*10" "${SCRIPT_DIR}/data/info_summarizer/non_numeric.out"; then
    echo "✗ Test failed: non_numeric - DP stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/non_numeric.out"
    exit 1
fi
if ! grep -q "^AF[[:space:]]*0.3333[[:space:]]*0.3[[:space:]]*0.5" "${SCRIPT_DIR}/data/info_summarizer/non_numeric.out"; then
    echo "✗ Test failed: non_numeric - AF stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/non_numeric.out"
    exit 1
fi

echo "✓ Test 6 passed"

###############################################################################
# Test 7: Flags => interpret as '1' if present
###############################################################################
echo "Test 7: Flags"
cat > "${SCRIPT_DIR}/data/info_summarizer/flags.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	SOMATIC;DP=10;FOO
1	200	rs2	C	T	.	PASS	SOMATIC;DP=20
1	300	.	G	A	.	PASS	AF=0.3
EOF

debugHex "${SCRIPT_DIR}/data/info_summarizer/flags.vcf"

"${SUMMARIZER_BIN}" --info "DP,SOMATIC,FOO" < "${SCRIPT_DIR}/data/info_summarizer/flags.vcf" \
  > "${SCRIPT_DIR}/data/info_summarizer/flags.out"

# DP => [10, 20] => mean=15, median=15, mode=10
# SOMATIC => [1, 1] => mean=1, median=1, mode=1
# FOO => [1]
if ! grep -q "^DP[[:space:]]*15.0000[[:space:]]*15[[:space:]]*10" "${SCRIPT_DIR}/data/info_summarizer/flags.out"; then
    echo "✗ Test failed: flags - DP mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/flags.out"
    exit 1
fi
if ! grep -q "^SOMATIC[[:space:]]*1.0000[[:space:]]*1[[:space:]]*1" "${SCRIPT_DIR}/data/info_summarizer/flags.out"; then
    echo "✗ Test failed: flags - SOMATIC mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/flags.out"
    exit 1
fi
if ! grep -q "^FOO[[:space:]]*1.0000[[:space:]]*1[[:space:]]*1" "${SCRIPT_DIR}/data/info_summarizer/flags.out"; then
    echo "✗ Test failed: flags - FOO mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/flags.out"
    exit 1
fi

echo "✓ Test 7 passed"

###############################################################################
# Test 8: Malformed line
###############################################################################
echo "Test 8: Malformed line"
cat > "${SCRIPT_DIR}/data/info_summarizer/invalid.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10
ThisIsInvalid
1	200	rs2	C	T	.	PASS	DP=20
EOF

debugHex "${SCRIPT_DIR}/data/info_summarizer/invalid.vcf"

"${SUMMARIZER_BIN}" --info "DP" \
  < "${SCRIPT_DIR}/data/info_summarizer/invalid.vcf" \
  2> "${SCRIPT_DIR}/data/info_summarizer/invalid.err" \
  > "${SCRIPT_DIR}/data/info_summarizer/invalid.out"

if ! grep -q "Warning: Skipping malformed VCF line" "${SCRIPT_DIR}/data/info_summarizer/invalid.err"; then
    echo "✗ Test failed: invalid lines - missing warning"
    exit 1
fi

# DP => [10,20] => mean=15, median=15, mode=10
if ! grep -q "^DP[[:space:]]*15.0000[[:space:]]*15[[:space:]]*10" "${SCRIPT_DIR}/data/info_summarizer/invalid.out"; then
    echo "✗ Test failed: invalid lines - DP stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/invalid.out"
    exit 1
fi

echo "✓ Test 8 passed"

###############################################################################
# Test 9: No #CHROM => error
###############################################################################
echo "Test 9: No #CHROM"
cat > "${SCRIPT_DIR}/data/info_summarizer/no_chrom.vcf" <<EOF
1	100	.	A	G	.	PASS	DP=10
1	200	rs2	C	T	.	PASS	DP=20
EOF

debugHex "${SCRIPT_DIR}/data/info_summarizer/no_chrom.vcf"

if "${SUMMARIZER_BIN}" --info "DP" \
   < "${SCRIPT_DIR}/data/info_summarizer/no_chrom.vcf" \
   2> "${SCRIPT_DIR}/data/info_summarizer/no_chrom.err" || true
then
    :
fi

if ! grep -q "Error: VCF header (#CHROM) not found before records" \
      "${SCRIPT_DIR}/data/info_summarizer/no_chrom.err"; then
    echo "✗ Test failed: no_chrom - expected error not found"
    cat "${SCRIPT_DIR}/data/info_summarizer/no_chrom.err"
    exit 1
fi

echo "✓ Test 9 passed"

###############################################################################
# Test 10: Partial final line (no trailing newline)
###############################################################################
echo "Test 10: Partial final line"
cat > "${SCRIPT_DIR}/data/info_summarizer/partial_unix.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=5
1	200	rs2	T	C	.	PASS	DP=15
EOF

truncate -s -1 "${SCRIPT_DIR}/data/info_summarizer/partial_unix.vcf"

debugHex "${SCRIPT_DIR}/data/info_summarizer/partial_unix.vcf"

"${SUMMARIZER_BIN}" --info "DP" \
  < "${SCRIPT_DIR}/data/info_summarizer/partial_unix.vcf" \
  > "${SCRIPT_DIR}/data/info_summarizer/partial_unix.out"

# DP => [5,15] => mean=10, median=10, mode=5
if ! grep -q "^DP[[:space:]]*10.0000[[:space:]]*10[[:space:]]*5" \
   "${SCRIPT_DIR}/data/info_summarizer/partial_unix.out"; then
    echo "✗ Test failed: partial line - DP stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/partial_unix.out"
    exit 1
fi

echo "✓ Test 10 passed"

###############################################################################
# Test 11: Windows CRLF line endings
###############################################################################
echo "Test 11: CRLF"
cat > "${SCRIPT_DIR}/data/info_summarizer/crlf_unix.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	.	PASS	DP=10
1	200	.	T	C	.	PASS	DP=30
EOF

sed 's/$/\r/' "${SCRIPT_DIR}/data/info_summarizer/crlf_unix.vcf" \
  > "${SCRIPT_DIR}/data/info_summarizer/crlf.vcf"

debugHex "${SCRIPT_DIR}/data/info_summarizer/crlf.vcf"

"${SUMMARIZER_BIN}" --info "DP" \
  < "${SCRIPT_DIR}/data/info_summarizer/crlf.vcf" \
  > "${SCRIPT_DIR}/data/info_summarizer/crlf.out"

# DP => [10,30] => mean=20, median=20, mode=10
if ! grep -q "^DP[[:space:]]*20.0000[[:space:]]*20[[:space:]]*10" \
   "${SCRIPT_DIR}/data/info_summarizer/crlf.out"
then
    echo "✗ Test failed: crlf - DP stats mismatch"
    cat "${SCRIPT_DIR}/data/info_summarizer/crlf.out"
    exit 1
fi

echo "✓ Test 11 passed"

###############################################################################
echo "✓ All tests passed"
