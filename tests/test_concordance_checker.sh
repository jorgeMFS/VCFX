#!/bin/bash

# Test script for VCFX_concordance_checker

# Stop on first error
set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# The root directory is one level up from the tests directory
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Create test data directory if it doesn't exist
mkdir -p data

# Create a temporary directory for test outputs
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Create test VCF file
cat > data/concordance_input.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	T	100	PASS	.	GT	0/1	0/1
1	200	rs2	G	C	100	PASS	.	GT	1/1	1/1
1	300	rs3	T	A	100	PASS	.	GT	0/0	0/1
1	400	rs4	C	G	100	PASS	.	GT	1/1	0/0
1	500	rs5	A	T,C	100	PASS	.	GT	1/2	1/2
1	600	rs6	G	C,A	100	PASS	.	GT	2/1	1/2
1	700	rs7	T	A,G	100	PASS	.	GT	./.	./.
1	800	rs8	C	G,T	100	PASS	.	GT	0/1	./. 
EOF

# Create expected output file
cat > data/concordance_expected.tsv << 'EOF'
CHROM	POS	ID	REF	ALT	SAMPLE1_GT	SAMPLE2_GT	Concordance
1	100	rs1	A	T	0/1	0/1	Concordant
1	200	rs2	G	C	1/1	1/1	Concordant
1	300	rs3	T	A	0/0	0/1	Discordant
1	400	rs4	C	G	1/1	0/0	Discordant
1	500	rs5	A	T,C	1/2	1/2	Concordant
1	600	rs6	G	C,A	1/2	1/2	Concordant
EOF

CONCORDANCE_CHECKER="../build/src/VCFX_concordance_checker/VCFX_concordance_checker"

# Test 1: Basic concordance checking
echo "Running test: Basic concordance checking"
"$CONCORDANCE_CHECKER" --samples "SAMPLE1 SAMPLE2" < data/concordance_input.vcf > "${TEMP_DIR}/output1.tsv"
if diff -w data/concordance_expected.tsv "${TEMP_DIR}/output1.tsv" > /dev/null; then
    echo "✓ Test passed"
else
    echo "✗ Test failed: Output differs from expected"
    echo "Differences:"
    diff -w data/concordance_expected.tsv "${TEMP_DIR}/output1.tsv"
    echo "Some tests for VCFX_concordance_checker failed"
    exit 1
fi

# Test 2: Multi-allelic variant handling
echo "Running test: Multi-allelic variant handling"
if grep -q "1	500	rs5	A	T,C" "${TEMP_DIR}/output1.tsv" && \
   grep -q "1	600	rs6	G	C,A" "${TEMP_DIR}/output1.tsv"; then
    echo "✓ Test passed"
else
    echo "✗ Test failed: Multi-allelic variants not handled correctly"
    echo "Some tests for VCFX_concordance_checker failed"
    exit 1
fi

# Test 3: Missing data handling
echo "Running test: Missing data handling"
if grep -q "1	700" "${TEMP_DIR}/output1.tsv"; then
    echo "✗ Test failed: Missing data not handled correctly"
    echo "Some tests for VCFX_concordance_checker failed"
    exit 1
else
    echo "✓ Test passed"
fi

# Test 4: Error handling for missing samples argument
echo "Running test: Error handling for missing samples argument"
set +e  # Disable exit on error for these tests
"$CONCORDANCE_CHECKER" < data/concordance_input.vcf > "${TEMP_DIR}/output4.tsv" 2>&1
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ] && grep -q "samples" "${TEMP_DIR}/output4.tsv"; then
    echo "✓ Test passed"
else
    echo "✗ Test failed: Missing samples error not handled correctly"
    echo "Some tests for VCFX_concordance_checker failed"
    exit 1
fi

# Test 5: Error handling for sample not found
echo "Running test: Error handling for sample not found"
"$CONCORDANCE_CHECKER" --samples "SAMPLE1 NONEXISTENT" < data/concordance_input.vcf > "${TEMP_DIR}/output5.tsv" 2>&1
EXIT_CODE=$?
set -e  # Re-enable exit on error
if [ $EXIT_CODE -ne 0 ] && grep -q "not found" "${TEMP_DIR}/output5.tsv"; then
    echo "✓ Test passed"
else
    echo "✗ Test failed: Sample not found error not handled correctly"
    echo "Some tests for VCFX_concordance_checker failed"
    exit 1
fi

# Test 6: Help message
echo "Running test: Help message"
set +e  # Disable exit on error for help test
"$CONCORDANCE_CHECKER" --help > "${TEMP_DIR}/help.txt" 2>&1
set -e  # Re-enable exit on error

# Check if the help message contains expected content
if grep -q "VCFX_concordance_checker" "${TEMP_DIR}/help.txt" && \
   grep -q "Usage:" "${TEMP_DIR}/help.txt" && \
   grep -q "Options:" "${TEMP_DIR}/help.txt"; then
    echo "✓ Test passed"
else
    echo "✗ Test failed: Help message does not contain expected content"
    echo "Help message content:"
    cat "${TEMP_DIR}/help.txt"
    echo "Some tests for VCFX_concordance_checker failed"
    exit 1
fi

echo "All tests for VCFX_concordance_checker passed"
exit 0 