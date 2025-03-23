#!/usr/bin/env bash

# Test script for VCFX_nonref_filter
# This script tests that the tool correctly filters out variants where all samples are homozygous reference

# Clear previous outputs
mkdir -p out/nonref_filter

# Colors for test results
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'  # No Color

# Set the executable
EXEC="../build/src/VCFX_nonref_filter/VCFX_nonref_filter"

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "❌ Error: $EXEC not found. Make sure to build the project first."
    exit 1
fi

# Function to run a test and compare output
run_test() {
    local test_name=$1
    local input_file="data/nonref_filter/${test_name}.vcf"
    local expected_file="expected/nonref_filter/${test_name}_out.vcf"
    local output_file="out/nonref_filter/${test_name}.out"
    
    echo "Running test: $test_name"
    $EXEC < "$input_file" > "$output_file"
    
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
    local output_file="out/nonref_filter/help.out"
    $EXEC --help > "$output_file"
    
    # Check if help output contains expected text
    if grep -q "VCFX_nonref_filter: Exclude variants if all samples are homozygous reference" "$output_file"; then
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
echo "=== Testing VCFX_nonref_filter ==="

# Initialize counters
tests_passed=0
tests_failed=0
total_tests=0

# Run tests
for test in "basic" "complex" "malformed" "no_gt"; do
    ((total_tests++))
    if run_test "$test"; then
        ((tests_passed++))
    else
        ((tests_failed++))
    fi
done

# Test help message
((total_tests++))
if test_help_message; then
    ((tests_passed++))
else
    ((tests_failed++))
fi

# Report results
echo "Tests completed: $tests_passed passed, $tests_failed failed (out of $total_tests total)"

if [ $tests_failed -eq 0 ]; then
    echo -e "${GREEN}✅ All VCFX_nonref_filter tests passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ Some VCFX_nonref_filter tests failed.${NC}"
    exit 1
fi 