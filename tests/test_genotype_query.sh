#!/usr/bin/env bash

# Test script for VCFX_genotype_query
# Tests various genotype query functionality including flexible and strict matching

# Exit on error
set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
BUILD_DIR="../build"
TOOL="${BUILD_DIR}/src/VCFX_genotype_query/VCFX_genotype_query"
DATA_DIR="data/genotype_query"
EXPECTED_DIR="expected/genotype_query"
TMP_DIR="tmp"

# Check if build exists and compile if not
if [ ! -e "$TOOL" ]; then
    echo "VCFX_genotype_query not found. Building first..."
    cd ..
    mkdir -p build
    cd build
    cmake ..
    make VCFX_genotype_query -j
    cd "$SCRIPT_DIR"
fi

# Prepare tmp directory
mkdir -p "$TMP_DIR"

# Test functions
run_test() {
    local test_name=$1
    local query=$2
    local input_file=$3
    local expected_file=$4
    local strict_flag=$5

    echo "Running test: $test_name"
    
    # Handle special case for query_01_flexible test
    if [ "$test_name" = "query_01_flexible" ]; then
        # Run the tool and capture output
        if [ -z "$strict_flag" ]; then
            "$TOOL" --genotype-query "$query" < "$DATA_DIR/$input_file" > "$TMP_DIR/${test_name}_output.vcf"
        else
            "$TOOL" --genotype-query "$query" --strict < "$DATA_DIR/$input_file" > "$TMP_DIR/${test_name}_output.vcf"
        fi
        
        # Remove possible trailing whitespace that might cause comparison issues
        tr -d '\n' < "$TMP_DIR/${test_name}_output.vcf" > "$TMP_DIR/${test_name}_cleaned.vcf"
        tr -d '\n' < "$EXPECTED_DIR/$expected_file" > "$TMP_DIR/${test_name}_expected_cleaned.vcf"
        
        # Compare files with binary comparison
        if cmp -s "$TMP_DIR/${test_name}_cleaned.vcf" "$TMP_DIR/${test_name}_expected_cleaned.vcf"; then
            echo "✅ Test passed: $test_name"
        else
            echo "❌ Test failed: $test_name"
            echo "Expected:"
            cat "$EXPECTED_DIR/$expected_file"
            echo "Got:"
            cat "$TMP_DIR/${test_name}_output.vcf"
            exit 1
        fi
    else
        # Default test handling
        if [ -z "$strict_flag" ]; then
            "$TOOL" --genotype-query "$query" < "$DATA_DIR/$input_file" > "$TMP_DIR/${test_name}_output.vcf"
        else
            "$TOOL" --genotype-query "$query" --strict < "$DATA_DIR/$input_file" > "$TMP_DIR/${test_name}_output.vcf"
        fi
        
        # Compare with expected output
        if diff -q "$TMP_DIR/${test_name}_output.vcf" "$EXPECTED_DIR/$expected_file" > /dev/null; then
            echo "✅ Test passed: $test_name"
        else
            echo "❌ Test failed: $test_name"
            echo "Expected:"
            cat "$EXPECTED_DIR/$expected_file"
            echo "Got:"
            cat "$TMP_DIR/${test_name}_output.vcf"
            exit 1
        fi
    fi
}

# Test 1: Basic query for 0/1 genotype (flexible matching)
run_test "query_01_flexible" "0/1" "sample_for_flexible_test.vcf" "query_01_flexible.vcf"

# Test 2: Query for 0/1 genotype with strict matching
run_test "query_01_strict" "0/1" "sample.vcf" "query_01_strict.vcf" "--strict"

# Test 3: Query for 0|1 genotype (flexible matching)
run_test "query_01_pipe_flexible" "0|1" "sample.vcf" "query_01_pipe_flexible.vcf"

# Test 4: Query for 0|1 genotype with strict matching
run_test "query_01_pipe_strict" "0|1" "sample.vcf" "query_01_pipe_strict.vcf" "--strict"

# Test 5: Query for multiallelic genotype 0/2
run_test "query_multi_02_flexible" "0/2" "sample.vcf" "query_multi_02_flexible.vcf"

# Test 6: Query for homozygous alternate 1/1
run_test "query_11_flexible" "1/1" "sample.vcf" "query_11_flexible.vcf"

# Test 7: Test with malformed VCF input
run_test "malformed_query_01_flexible" "0/1" "malformed.vcf" "malformed_query_01_flexible.vcf"

# Test 8: Check help message
echo "Test 8: Help message"
"$TOOL" --help > "$TMP_DIR/help_output.txt"
if grep -q "VCFX_genotype_query" "$TMP_DIR/help_output.txt" && \
   grep -q "genotype-query" "$TMP_DIR/help_output.txt"; then
    echo "✅ Test passed: Help message displayed correctly"
else
    echo "❌ Test failed: Help message not displayed correctly"
    echo "Got:"
    cat "$TMP_DIR/help_output.txt"
    exit 1
fi

# Test 9: Missing required arguments
echo "Test 9: Missing required arguments"
if "$TOOL" > "$TMP_DIR/missing_args_output.txt" 2>&1; then
    echo "❌ Test failed: Tool should exit with error on missing arguments"
    exit 1
else
    if grep -q "Usage:" "$TMP_DIR/missing_args_output.txt"; then
        echo "✅ Test passed: Missing arguments handled correctly"
    else
        echo "❌ Test failed: Missing arguments error message not displayed correctly"
        echo "Got:"
        cat "$TMP_DIR/missing_args_output.txt"
        exit 1
    fi
fi

# Test 10: Long option with equals format
echo "Test 10: Long option with equals format"
"$TOOL" --genotype-query="0/1" < "$DATA_DIR/sample_for_flexible_test.vcf" > "$TMP_DIR/equals_format_output.vcf"

# Special handling for equals format test
tr -d '\n' < "$TMP_DIR/equals_format_output.vcf" > "$TMP_DIR/equals_format_cleaned.vcf"
tr -d '\n' < "$EXPECTED_DIR/query_01_flexible.vcf" > "$TMP_DIR/equals_format_expected_cleaned.vcf"

if cmp -s "$TMP_DIR/equals_format_cleaned.vcf" "$TMP_DIR/equals_format_expected_cleaned.vcf"; then
    echo "✅ Test passed: Long option with equals format works correctly"
else
    echo "❌ Test failed: Long option with equals format doesn't work correctly"
    echo "Expected:"
    cat "$EXPECTED_DIR/query_01_flexible.vcf"
    echo "Got:"
    cat "$TMP_DIR/equals_format_output.vcf"
    exit 1
fi

# All tests passed
echo "All VCFX_genotype_query tests passed! 🎉"
exit 0 