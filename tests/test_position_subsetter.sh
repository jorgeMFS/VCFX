#!/usr/bin/env bash

# Test script for VCFX_position_subsetter
# This script tests various aspects of the position subsetter functionality

# Set up variables
EXEC="../build/src/VCFX_position_subsetter/VCFX_position_subsetter"
failures=0

# Check if the executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: VCFX_position_subsetter executable not found at $EXEC"
    echo "Please build the project before running this test."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p out/position_subsetter

echo "=== Testing VCFX_position_subsetter ==="

# Function to run test and check result with diff
run_test() {
    test_num=$1
    test_name=$2
    vcf_file=$3
    region=$4
    expected_contents=$5
    
    echo "Test $test_num: $test_name"
    $EXEC --region "$region" < "$vcf_file" > "out/position_subsetter/test${test_num}.vcf" 2> "out/position_subsetter/test${test_num}.err"
    exit_code=$?
    
    if [ $exit_code -ne 0 ] && [ "$6" != "expect_failure" ]; then
        echo "  Test $test_num failed: command exited with code $exit_code"
        cat "out/position_subsetter/test${test_num}.err"
        ((failures++))
    elif [ $exit_code -eq 0 ] && [ "$6" == "expect_failure" ]; then
        echo "  Test $test_num failed: command succeeded but failure was expected"
        ((failures++))
    elif [ "$6" != "expect_failure" ]; then
        # Use grep to check if all expected lines are in the output
        missing_lines=0
        while IFS= read -r line; do
            if ! grep -q -F "$line" "out/position_subsetter/test${test_num}.vcf"; then
                echo "  Missing line in output: $line"
                ((missing_lines++))
            fi
        done <<< "$expected_contents"
        
        # Check if there are unexpected extra lines in the output
        output_lines=$(grep -v -e '^$' "out/position_subsetter/test${test_num}.vcf" | wc -l)
        expected_lines=$(echo "$expected_contents" | grep -v -e '^$' | wc -l)
        
        if [ $missing_lines -gt 0 ] || [ $output_lines -ne $expected_lines ]; then
            echo "  Test $test_num failed: output does not match expected"
            echo "  Expected:"
            echo "$expected_contents"
            echo "  Actual:"
            cat "out/position_subsetter/test${test_num}.vcf"
            ((failures++))
        else
            echo "  Test $test_num passed."
        fi
    else
        echo "  Test $test_num passed (expected failure)."
    fi
}

# Function to test help display
test_help() {
    test_num=$1
    echo "Test $test_num: Help display"
    
    $EXEC --help > "out/position_subsetter/test${test_num}.out" 2>&1
    if ! grep -q "Usage:" "out/position_subsetter/test${test_num}.out"; then
        echo "  Test $test_num failed: help message does not contain 'Usage:'"
        ((failures++))
    elif ! grep -q "VCFX_position_subsetter" "out/position_subsetter/test${test_num}.out"; then
        echo "  Test $test_num failed: help message does not contain tool name"
        ((failures++))
    else
        echo "  Test $test_num passed."
    fi
}

# Test 1: Basic functionality - chr1:100-200
expected1="##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	50	PASS	DP=10
chr1	150	.	C	G	60	PASS	DP=15
chr1	200	.	G	A	70	PASS	DP=20"
run_test 1 "Basic subsetting chr1:100-200" \
    "data/position_subsetter/multi_chrom.vcf" \
    "chr1:100-200" \
    "$expected1"

# Test 2: Partial match - chr2:150-250
expected2="##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr2	200	.	G	T	65	PASS	DP=18"
run_test 2 "Partial match subsetting chr2:150-250" \
    "data/position_subsetter/multi_chrom.vcf" \
    "chr2:150-250" \
    "$expected2"

# Test 3: Range at the end - chr3:250-400
expected3="##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr3	250	.	T	G	55	PASS	DP=11
chr3	350	.	G	C	65	PASS	DP=13"
run_test 3 "Range at the end - chr3:250-400" \
    "data/position_subsetter/multi_chrom.vcf" \
    "chr3:250-400" \
    "$expected3"

# Test 4: No matching chromosome - chr4:100-200
expected4="##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
run_test 4 "No matching chromosome - chr4:100-200" \
    "data/position_subsetter/multi_chrom.vcf" \
    "chr4:100-200" \
    "$expected4"

# Test 5: Invalid region format
run_test 5 "Invalid region format" \
    "data/position_subsetter/multi_chrom.vcf" \
    "chr1:100" \
    "$expected4" \
    "expect_failure"

# Test 6: End smaller than start
run_test 6 "End smaller than start" \
    "data/position_subsetter/multi_chrom.vcf" \
    "chr1:200-100" \
    "$expected4" \
    "expect_failure"

# Test 7: Malformed input VCF
expected7="##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	150	.	C	G	60	PASS	DP=15"
run_test 7 "Malformed input VCF" \
    "data/position_subsetter/malformed.vcf" \
    "chr1:100-200" \
    "$expected7"

# Test 8: Missing header line
expected8="##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
run_test 8 "Missing header line" \
    "data/position_subsetter/no_header.vcf" \
    "chr1:100-200" \
    "$expected8"

# Test 9: Help message
test_help 9

# Report results
if [ $failures -eq 0 ]; then
    echo "All VCFX_position_subsetter tests passed!"
    exit 0
else
    echo "$failures VCFX_position_subsetter tests failed."
    exit 1
fi 