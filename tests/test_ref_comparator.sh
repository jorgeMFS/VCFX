#!/usr/bin/env bash

# Test script for VCFX_ref_comparator
# This script tests various aspects of the reference comparator functionality

# Set up variables
EXEC="../build/src/VCFX_ref_comparator/VCFX_ref_comparator"
REF_FASTA="data/ref_comparator/test_reference.fasta"
failures=0

# Check if the executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: VCFX_ref_comparator executable not found at $EXEC"
    echo "Please build the project before running this test."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p out/ref_comparator

echo "=== Testing VCFX_ref_comparator ==="

# Function to run test and check result
run_test() {
    test_num=$1
    test_name=$2
    vcf_file=$3
    expected_file=$4
    
    echo "Test $test_num: $test_name"
    $EXEC --reference "$REF_FASTA" < "$vcf_file" > "out/ref_comparator/test${test_num}.vcf" 2> "out/ref_comparator/test${test_num}.err"
    exit_code=$?
    
    if [ $exit_code -ne 0 ]; then
        echo "  Test $test_num failed: command exited with code $exit_code"
        cat "out/ref_comparator/test${test_num}.err"
        ((failures++))
    else
        if ! diff -q "$expected_file" "out/ref_comparator/test${test_num}.vcf" > /dev/null; then
            echo "  Test $test_num failed: output does not match expected"
            echo "  Expected:"
            cat "$expected_file"
            echo "  Actual:"
            cat "out/ref_comparator/test${test_num}.vcf"
            ((failures++))
        else
            echo "  Test $test_num passed."
        fi
    fi
}

# Function to test help display
test_help() {
    test_num=$1
    echo "Test $test_num: Help display"
    
    $EXEC --help > "out/ref_comparator/test${test_num}.out" 2>&1
    if ! grep -q "Usage:" "out/ref_comparator/test${test_num}.out"; then
        echo "  Test $test_num failed: help message does not contain 'Usage:'"
        ((failures++))
    elif ! grep -q "VCFX_ref_comparator" "out/ref_comparator/test${test_num}.out"; then
        echo "  Test $test_num failed: help message does not contain tool name"
        ((failures++))
    else
        echo "  Test $test_num passed."
    fi
}

# Function to test invalid reference file
test_invalid_ref() {
    test_num=$1
    echo "Test $test_num: Invalid reference file"
    
    $EXEC --reference "non_existent_file.fasta" < "data/ref_comparator/basic.vcf" > "out/ref_comparator/test${test_num}.vcf" 2> "out/ref_comparator/test${test_num}.err"
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "  Test $test_num failed: command succeeded with non-existent reference file"
        ((failures++))
    else
        if ! grep -q "Error: cannot open reference" "out/ref_comparator/test${test_num}.err"; then
            echo "  Test $test_num failed: error message does not contain expected text"
            ((failures++))
        else
            echo "  Test $test_num passed (expected failure with proper error message)."
        fi
    fi
}

# Test 1: Basic functionality - valid CHROM, POS, REF, ALT
run_test 1 "Basic functionality" \
    "data/ref_comparator/basic.vcf" \
    "expected/ref_comparator/basic_out.vcf"

# Test 2: ALT matches REF in some cases
run_test 2 "ALT matches REF" \
    "data/ref_comparator/alt_matches_ref.vcf" \
    "expected/ref_comparator/alt_matches_ref_out.vcf"

# Test 3: Invalid positions and unknown chromosomes
run_test 3 "Invalid positions and chromosomes" \
    "data/ref_comparator/invalid_positions.vcf" \
    "expected/ref_comparator/invalid_positions_out.vcf"

# Test 4: Malformed input
run_test 4 "Malformed input" \
    "data/ref_comparator/malformed.vcf" \
    "expected/ref_comparator/malformed_out.vcf"

# Test 5: Help message
test_help 5

# Test 6: Invalid reference file
test_invalid_ref 6

# Performance test with large input would be added here

# Report results
if [ $failures -eq 0 ]; then
    echo "All VCFX_ref_comparator tests passed!"
    exit 0
else
    echo "$failures VCFX_ref_comparator tests failed."
    exit 1
fi 