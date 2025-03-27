#!/bin/bash

# Test script for VCFX_compressor

# Stop on first error
set -e

# Change to the directory containing this script
cd "$(dirname "$0")"

# Get the root directory
ROOT_DIR="$(cd .. && pwd)"

# Make sure the build directory exists and the tool is built
if [ ! -f "${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor" ]; then
    echo "Building VCFX_compressor..."
    mkdir -p "${ROOT_DIR}/build"
    cd "${ROOT_DIR}/build"
    cmake ..
    make VCFX_compressor
    cd "${ROOT_DIR}/tests"
fi

# Create a temporary directory for test outputs
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Create test data directory if it doesn't exist
mkdir -p data

# Create test VCF file
cat > data/compressor_input.vcf << 'EOF'
##fileformat=VCFv4.2
##source=VCFX_test
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs1	A	G	100	PASS	DP=50;AF=0.1
chr1	200	rs2	C	T	80	PASS	DP=60;AF=0.2
chr1	300	rs3	G	A	60	PASS	DP=70;AF=0.3
chr2	400	rs4	T	C	90	PASS	DP=80;AF=0.4
chr2	500	rs5	A	G	70	PASS	DP=90;AF=0.5
EOF

# Create a large VCF file to test buffer handling
cat > data/compressor_large.vcf << 'EOF'
##fileformat=VCFv4.2
##source=VCFX_test
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

# Add 1000 lines to test buffer handling
for i in {1..1000}; do
    echo "chr1	$i	rs$i	A	G	100	PASS	DP=$i" >> data/compressor_large.vcf
done

# Create an empty VCF file
touch data/compressor_empty.vcf

# Test 1: Basic compression and decompression
echo "Test 1: Basic compression and decompression"
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --compress < data/compressor_input.vcf > ${TEMP_DIR}/compressed.vcf.gz
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --decompress < ${TEMP_DIR}/compressed.vcf.gz > ${TEMP_DIR}/decompressed.vcf

if diff -q data/compressor_input.vcf ${TEMP_DIR}/decompressed.vcf > /dev/null; then
    echo "  ✓ Basic compression/decompression test passed"
else
    echo "  ✗ Basic compression/decompression test failed"
    echo "Differences:"
    diff data/compressor_input.vcf ${TEMP_DIR}/decompressed.vcf
    exit 1
fi

# Test 2: Large file compression and decompression
echo "Test 2: Large file compression and decompression"
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --compress < data/compressor_large.vcf > ${TEMP_DIR}/compressed_large.vcf.gz
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --decompress < ${TEMP_DIR}/compressed_large.vcf.gz > ${TEMP_DIR}/decompressed_large.vcf

if diff -q data/compressor_large.vcf ${TEMP_DIR}/decompressed_large.vcf > /dev/null; then
    echo "  ✓ Large file compression/decompression test passed"
else
    echo "  ✗ Large file compression/decompression test failed"
    echo "Differences:"
    diff data/compressor_large.vcf ${TEMP_DIR}/decompressed_large.vcf
    exit 1
fi

# Test 3: Empty file handling
echo "Test 3: Empty file handling"
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --compress < data/compressor_empty.vcf > ${TEMP_DIR}/compressed_empty.vcf.gz
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --decompress < ${TEMP_DIR}/compressed_empty.vcf.gz > ${TEMP_DIR}/decompressed_empty.vcf

if diff -q data/compressor_empty.vcf ${TEMP_DIR}/decompressed_empty.vcf > /dev/null; then
    echo "  ✓ Empty file handling test passed"
else
    echo "  ✗ Empty file handling test failed"
    echo "Differences:"
    diff data/compressor_empty.vcf ${TEMP_DIR}/decompressed_empty.vcf
    exit 1
fi

# Test 4: Error handling - invalid compressed data
echo "Test 4: Error handling - invalid compressed data"
set +e  # Don't exit on error for this test
echo "This is not a valid gzip file" > ${TEMP_DIR}/invalid.gz
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --decompress < ${TEMP_DIR}/invalid.gz > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: inflate failed" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Invalid compressed data error handling test passed"
else
    echo "  ✗ Invalid compressed data error handling test failed"
    echo "Expected non-zero exit code and error message about invalid data"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 5: Error handling - no mode specified
echo "Test 5: Error handling - no mode specified"
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor < data/compressor_input.vcf > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: must specify exactly one of --compress or --decompress" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ No mode specified error handling test passed"
else
    echo "  ✗ No mode specified error handling test failed"
    echo "Expected non-zero exit code and error message about mode specification"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 6: Error handling - both modes specified
echo "Test 6: Error handling - both modes specified"
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --compress --decompress < data/compressor_input.vcf > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: must specify exactly one of --compress or --decompress" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Both modes specified error handling test passed"
else
    echo "  ✗ Both modes specified error handling test failed"
    echo "Expected non-zero exit code and error message about mode specification"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 7: Help message
echo "Test 7: Help message"
${ROOT_DIR}/build/src/VCFX_compressor/VCFX_compressor --help > ${TEMP_DIR}/help_output.txt

if grep -q "VCFX_compressor" ${TEMP_DIR}/help_output.txt && \
   grep -q "Usage:" ${TEMP_DIR}/help_output.txt && \
   grep -q "Options:" ${TEMP_DIR}/help_output.txt; then
    echo "  ✓ Help message test passed"
else
    echo "  ✗ Help message test failed"
    echo "Expected help message containing tool description"
    echo "Got output:"
    cat ${TEMP_DIR}/help_output.txt
    exit 1
fi

echo "All tests for VCFX_compressor passed!" 