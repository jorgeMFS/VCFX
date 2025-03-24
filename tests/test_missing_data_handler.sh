#!/bin/bash

# Test script for VCFX_missing_data_handler

# Stop on first error
set -e

# Change to the directory containing this script
cd "$(dirname "$0")"

# Get the root directory
ROOT_DIR="$(cd .. && pwd)"

# Make sure the build directory exists and the tool is built
if [ ! -f "${ROOT_DIR}/build/src/VCFX_missing_data_handler/VCFX_missing_data_handler" ]; then
    echo "Building VCFX_missing_data_handler..."
    mkdir -p "${ROOT_DIR}/build"
    cd "${ROOT_DIR}/build"
    cmake ..
    make VCFX_missing_data_handler
    cd "${ROOT_DIR}/tests"
fi

# Create a temporary directory for test outputs
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Test scenario 1: Default behavior (flag only)
echo "Test 1: Default behavior (flag only)"
${ROOT_DIR}/build/src/VCFX_missing_data_handler/VCFX_missing_data_handler < data/missing_data_handler_input.vcf > ${TEMP_DIR}/scenario1_output.vcf

# Compare with expected output
if diff -w -q ${TEMP_DIR}/scenario1_output.vcf expected/missing_data_handler_flagged.vcf > /dev/null; then
    echo "  ✓ Default behavior test passed"
else
    echo "  ✗ Default behavior test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/scenario1_output.vcf expected/missing_data_handler_flagged.vcf
    exit 1
fi

# Test scenario 2: Fill missing with default value (./)
echo "Test 2: Fill missing with default value (./.) "
${ROOT_DIR}/build/src/VCFX_missing_data_handler/VCFX_missing_data_handler --fill-missing < data/missing_data_handler_input.vcf > ${TEMP_DIR}/scenario2_output.vcf

# Compare with expected output
if diff -w -q ${TEMP_DIR}/scenario2_output.vcf expected/missing_data_handler_filled_default.vcf > /dev/null; then
    echo "  ✓ Fill with default value test passed"
else
    echo "  ✗ Fill with default value test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/scenario2_output.vcf expected/missing_data_handler_filled_default.vcf
    exit 1
fi

# Test scenario 3: Fill missing with custom value (0/0)
echo "Test 3: Fill missing with custom value (0/0)"
${ROOT_DIR}/build/src/VCFX_missing_data_handler/VCFX_missing_data_handler --fill-missing --default-genotype "0/0" < data/missing_data_handler_input.vcf > ${TEMP_DIR}/scenario3_output.vcf

# Compare with expected output
if diff -w -q ${TEMP_DIR}/scenario3_output.vcf expected/missing_data_handler_filled_custom.vcf > /dev/null; then
    echo "  ✓ Fill with custom value test passed"
else
    echo "  ✗ Fill with custom value test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/scenario3_output.vcf expected/missing_data_handler_filled_custom.vcf
    exit 1
fi

# Test scenario 4: Error handling - non-existent file
echo "Test 4: Error handling - non-existent file"
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_missing_data_handler/VCFX_missing_data_handler non_existent_file.vcf > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "Error: cannot open file" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Error handling test passed"
else
    echo "  ✗ Error handling test failed"
    echo "Expected non-zero exit code and error message about non-existent file"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test scenario 5: Multiple input files
echo "Test 5: Multiple input files processing"
# Use the same input file twice to simulate multiple files
cp data/missing_data_handler_input.vcf ${TEMP_DIR}/file1.vcf
cp data/missing_data_handler_input.vcf ${TEMP_DIR}/file2.vcf

# Run with the two copies as input
${ROOT_DIR}/build/src/VCFX_missing_data_handler/VCFX_missing_data_handler --fill-missing --default-genotype "0/0" ${TEMP_DIR}/file1.vcf ${TEMP_DIR}/file2.vcf > ${TEMP_DIR}/multiple_files_output.vcf

# Check that the output has data from both files (should have more lines than a single file)
INPUT_LINES=$(wc -l < data/missing_data_handler_input.vcf)
OUTPUT_LINES=$(wc -l < ${TEMP_DIR}/multiple_files_output.vcf)

if [ $OUTPUT_LINES -gt $INPUT_LINES ]; then
    # Further check: make sure missing values were filled with 0/0
    if grep -q "0/0" ${TEMP_DIR}/multiple_files_output.vcf; then
        echo "  ✓ Multiple input files test passed"
    else
        echo "  ✗ Multiple input files test failed - missing values were not filled"
        exit 1
    fi
else
    echo "  ✗ Multiple input files test failed - output doesn't contain data from both files"
    echo "    Input file lines: $INPUT_LINES"
    echo "    Output file lines: $OUTPUT_LINES"
    exit 1
fi

echo "All tests for VCFX_missing_data_handler passed!" 