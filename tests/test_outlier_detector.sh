#!/usr/bin/env bash

# Test script for VCFX_outlier_detector
# This script tests the tool's ability to detect outliers in variants and samples

# Clear previous outputs
mkdir -p out/outlier_detector

# Colors for test results
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'  # No Color

# Set the executable
EXEC="../build/src/VCFX_outlier_detector/VCFX_outlier_detector"

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "❌ Error: $EXEC not found. Make sure to build the project first."
    exit 1
fi

# Function to run a test and compare output
run_test() {
    local test_name=$1
    local input_file=$2
    local expected_file=$3
    local output_file=$4
    local cmd_args=$5
    
    echo "Running test: $test_name"
    $EXEC $cmd_args < "$input_file" > "$output_file"
    
    # Compare output with expected
    if diff -b "$output_file" "$expected_file" > /dev/null; then
        echo -e "${GREEN}✅ Test passed: $test_name${NC}"
        return 0
    else
        echo -e "${RED}❌ Test failed: $test_name${NC}"
        echo "Expected output (from $expected_file):"
        cat "$expected_file"
        echo "Actual output (from $output_file):"
        cat "$output_file"
        echo "Diff:"
        diff -b "$output_file" "$expected_file" | grep -v "^---" | grep -v "^+++"
        return 1
    fi
}

# Test help message
test_help_message() {
    echo "Running test: help_message"
    local output_file="out/outlier_detector/help.out"
    $EXEC --help > "$output_file"
    
    # Check if help output contains expected text
    if grep -q "VCFX_outlier_detector: Identify outliers among variants or samples based on a numeric metric" "$output_file"; then
        echo -e "${GREEN}✅ Test passed: help_message${NC}"
        return 0
    else
        echo -e "${RED}❌ Test failed: help_message${NC}"
        echo "Help message does not contain expected text:"
        cat "$output_file"
        return 1
    fi
}

# Run all tests
echo "=== Testing VCFX_outlier_detector ==="

# Initialize counters
tests_passed=0
tests_failed=0
total_tests=0

# Test 1: Variant mode with AF threshold
((total_tests++))
if run_test "variant_mode_af" \
    "data/outlier_detector/variant_mode.vcf" \
    "expected/outlier_detector/variant_af_0.1.txt" \
    "out/outlier_detector/variant_af_0.1.out" \
    "--metric AF --threshold 0.1 --variant"; then
    ((tests_passed++))
else
    ((tests_failed++))
fi

# Test 2: Variant mode with DP threshold
((total_tests++))
if run_test "variant_mode_dp" \
    "data/outlier_detector/variant_mode.vcf" \
    "expected/outlier_detector/variant_dp_150.txt" \
    "out/outlier_detector/variant_dp_150.out" \
    "--metric DP --threshold 150 --variant"; then
    ((tests_passed++))
else
    ((tests_failed++))
fi

# Test 3: Sample mode with GQ threshold
((total_tests++))
if run_test "sample_mode_gq" \
    "data/outlier_detector/sample_mode.vcf" \
    "expected/outlier_detector/sample_gq_35.txt" \
    "out/outlier_detector/sample_gq_35.out" \
    "--metric GQ --threshold 35 --sample"; then
    ((tests_passed++))
else
    ((tests_failed++))
fi

# Test 4: Malformed input handling
((total_tests++))
if run_test "malformed_handling" \
    "data/outlier_detector/malformed.vcf" \
    "expected/outlier_detector/malformed_af_0.1.txt" \
    "out/outlier_detector/malformed_af_0.1.out" \
    "--metric AF --threshold 0.1 --variant"; then
    ((tests_passed++))
else
    ((tests_failed++))
fi

# Test 5: Help message
((total_tests++))
if test_help_message; then
    ((tests_passed++))
else
    ((tests_failed++))
fi

# Report results
echo "Tests completed: $tests_passed passed, $tests_failed failed (out of $total_tests total)"

if [ $tests_failed -eq 0 ]; then
    echo -e "${GREEN}✅ All VCFX_outlier_detector tests passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ Some VCFX_outlier_detector tests failed.${NC}"
    exit 1
fi 