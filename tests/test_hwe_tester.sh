#!/bin/bash

# Test script for VCFX_hwe_tester

# Stop on first error
set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# The root directory is one level up from the tests directory
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Create test data directory if it doesn't exist
mkdir -p "${SCRIPT_DIR}/data/hwe_tester"

# Function to run a test and check its output
run_test() {
    local test_name="$1"
    local input="$2"
    local expected="$3"
    local description="$4"
    
    echo "Running test: $test_name"
    echo "Description: $description"
    
    # Write input to a file
    echo -e "$input" > "${SCRIPT_DIR}/data/hwe_tester/${test_name}.vcf"
    
    # Run the test
    "${ROOT_DIR}/build/src/VCFX_hwe_tester/VCFX_hwe_tester" < "${SCRIPT_DIR}/data/hwe_tester/${test_name}.vcf" > "${SCRIPT_DIR}/data/hwe_tester/${test_name}.out"
    
    # Compare output with expected
    if diff "${SCRIPT_DIR}/data/hwe_tester/${test_name}.out" <(echo -e "$expected") > /dev/null; then
        echo "✓ Test passed: $test_name"
    else
        echo "✗ Test failed: $test_name"
        echo "Expected:"
        echo -e "$expected"
        echo "Got:"
        cat "${SCRIPT_DIR}/data/hwe_tester/${test_name}.out"
        exit 1
    fi
    echo
}

# Test 1: Basic HWE test with perfect equilibrium
echo "Test 1: Basic HWE test with perfect equilibrium"
input1="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3
1\t100\trs1\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1"
expected1="CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue
1\t100\trs1\tA\tG\t1.000000"
run_test "basic_hwe" "$input1" "$expected1" "Testing basic HWE calculation with perfect equilibrium"

# Test 2: Multi-allelic variant (should be skipped)
echo "Test 2: Multi-allelic variant"
input2="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3
1\t100\trs1\tA\tG,T\t100\tPASS\t.\tGT\t0/0\t0/1\t1/2"
expected2="CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue"
run_test "multi_allelic" "$input2" "$expected2" "Testing that multi-allelic variants are skipped"

# Test 3: Missing genotypes
echo "Test 3: Missing genotypes"
input3="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3
1\t100\trs1\tA\tG\t100\tPASS\t.\tGT\t./.\t0/1\t1/1"
expected3="CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue
1\t100\trs1\tA\tG\t1.000000"
run_test "missing_genotypes" "$input3" "$expected3" "Testing handling of missing genotypes"

# Test 4: Different genotype formats
echo "Test 4: Different genotype formats"
input4="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3
1\t100\trs1\tA\tG\t100\tPASS\t.\tGT\t0|0\t0/1\t1|1"
expected4="CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue
1\t100\trs1\tA\tG\t1.000000"
run_test "genotype_formats" "$input4" "$expected4" "Testing different genotype formats (0/0 vs 0|0)"

# Test 5: No GT field in FORMAT
echo "Test 5: No GT field in FORMAT"
input5="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3
1\t100\trs1\tA\tG\t100\tPASS\t.\tDP\t10\t20\t30"
expected5="CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue"
run_test "no_gt_field" "$input5" "$expected5" "Testing handling of missing GT field"

# Test 6: Help message
echo "Test 6: Help message"
help_output=$("${ROOT_DIR}/build/src/VCFX_hwe_tester/VCFX_hwe_tester" --help)
if echo "$help_output" | grep -q "VCFX_hwe_tester: Perform Hardy-Weinberg Equilibrium"; then
    echo "✓ Test passed: Help message"
else
    echo "✗ Test failed: Help message"
    exit 1
fi

# Test 7: Performance with large file
echo "Test 7: Performance with large file"
# Create a large test file with 1000 variants
{
    echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3"
    for i in {1..1000}; do
        echo "1\t$i\trs$i\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1"
    done
} > "${SCRIPT_DIR}/data/hwe_tester/large.vcf"

time "${ROOT_DIR}/build/src/VCFX_hwe_tester/VCFX_hwe_tester" < "${SCRIPT_DIR}/data/hwe_tester/large.vcf" > /dev/null

echo "All VCFX_hwe_tester tests passed! 🎉" 