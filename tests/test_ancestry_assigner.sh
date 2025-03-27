#!/bin/bash

# Test script for VCFX_ancestry_assigner

# Stop on first error
set -e

# Change to the directory containing this script
cd "$(dirname "$0")"

# Get the root directory
ROOT_DIR="$(cd .. && pwd)"

# Make sure the build directory exists and the tool is built
if [ ! -f "${ROOT_DIR}/build/src/VCFX_ancestry_assigner/VCFX_ancestry_assigner" ]; then
    echo "Building VCFX_ancestry_assigner..."
    mkdir -p "${ROOT_DIR}/build"
    cd "${ROOT_DIR}/build"
    cmake ..
    make VCFX_ancestry_assigner
    cd "${ROOT_DIR}/tests"
fi

# Create a temporary directory for test outputs
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Create test data directory if it doesn't exist
mkdir -p data

# Create test frequency file
cat > data/ancestry_freq.tsv << 'EOF'
CHROM	POS	REF	ALT	EUR	ASN	AFR
chr1	10000	A	G	0.1	0.2	0.3
chr1	20000	C	T	0.2	0.3	0.4
chr1	30000	G	A	0.3	0.4	0.5
chr2	15000	T	C	0.4	0.5	0.6
chr2	25000	G	T	0.5	0.6	0.7
chr3	5000	A	C	0.6	0.7	0.8
chr3	12500	G	A	0.7	0.8	0.9
EOF

# Create test VCF file
cat > data/ancestry_input.vcf << 'EOF'
##fileformat=VCFv4.2
##source=VCFX_test
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0/1:30:99	0/0:25:95	1/1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0/0:28:99	0/1:32:98	0/0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1/1:15:45	1/1:14:52	0/1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0/1:45:99	0/0:42:99	1/1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP	1/1:10	0/1:8	0/0:9
chr3	5000	rs103	A	C	95	PASS	DP=72	GT:DP:GQ	0/0:30:99	1/1:32:98	0/1:27:97
chr3	12500	rs104	G	A	88	PASS	DP=65	GT	0/0	1/1	0/1
EOF

# Create expected output
# Note: These ancestry assignments are based on a proper likelihood calculation
# that considers the probability of each genotype given the population frequencies.
# For each sample, we compute P(genotype|freq) = (1-f)^2 for 0/0, 2f(1-f) for 0/1,
# and f^2 for 1/1, then sum the log-likelihoods across all variants.
# This is more statistically sound than just looking at frequency ordering.
cat > expected/ancestry_assigner_output.txt << 'EOF'
SAMPLE1	EUR
SAMPLE2	AFR
SAMPLE3	ASN
EOF

# Test 1: Basic functionality
echo "Test 1: Basic ancestry assignment"
${ROOT_DIR}/build/src/VCFX_ancestry_assigner/VCFX_ancestry_assigner --assign-ancestry data/ancestry_freq.tsv < data/ancestry_input.vcf > ${TEMP_DIR}/output.txt

# Compare with expected output
if diff -w -q ${TEMP_DIR}/output.txt expected/ancestry_assigner_output.txt > /dev/null; then
    echo "  ✓ Basic ancestry assignment test passed"
else
    echo "  ✗ Basic ancestry assignment test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/output.txt expected/ancestry_assigner_output.txt
    exit 1
fi

# Test 2: Error handling - missing frequency file
echo "Test 2: Error handling - missing frequency file"
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_ancestry_assigner/VCFX_ancestry_assigner --assign-ancestry non_existent_file.tsv > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: Unable to open frequency file" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Error handling test passed"
else
    echo "  ✗ Error handling test failed"
    echo "Expected non-zero exit code and error message about missing file"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 3: Error handling - empty frequency file
echo "Test 3: Error handling - empty frequency file"
touch ${TEMP_DIR}/empty.tsv
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_ancestry_assigner/VCFX_ancestry_assigner --assign-ancestry ${TEMP_DIR}/empty.tsv > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: Frequency file is empty" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Empty file error handling test passed"
else
    echo "  ✗ Empty file error handling test failed"
    echo "Expected non-zero exit code and error message about empty file"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 4: Error handling - invalid frequency file format
echo "Test 4: Error handling - invalid frequency file format"
cat > ${TEMP_DIR}/invalid.tsv << 'EOF'
CHROM	POS	REF	ALT
EOF
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_ancestry_assigner/VCFX_ancestry_assigner --assign-ancestry ${TEMP_DIR}/invalid.tsv > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: Frequency header must have at least 5 columns" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Invalid format error handling test passed"
else
    echo "  ✗ Invalid format error handling test failed"
    echo "Expected non-zero exit code and error message about invalid format"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 5: Help message
echo "Test 5: Help message"
${ROOT_DIR}/build/src/VCFX_ancestry_assigner/VCFX_ancestry_assigner --help > ${TEMP_DIR}/help_output.txt

if grep -q "VCFX_ancestry_assigner: Assign samples to ancestral populations" ${TEMP_DIR}/help_output.txt; then
    echo "  ✓ Help message test passed"
else
    echo "  ✗ Help message test failed"
    echo "Expected help message containing tool description"
    echo "Got output:"
    cat ${TEMP_DIR}/help_output.txt
    exit 1
fi

echo "All tests for VCFX_ancestry_assigner passed!" 