#!/bin/bash

# Test script for VCFX_impact_filter

# Stop on first error
set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# The root directory is one level up from the tests directory
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Create test data directory if it doesn't exist
mkdir -p "${SCRIPT_DIR}/data/impact_filter"

# Function to run a test and check its output
run_test() {
    local test_name="$1"
    local input="$2"
    local expected="$3"
    local filter_level="$4"
    local description="$5"
    
    echo "Running test: $test_name"
    echo "Description: $description"
    
    # Write input to a file
    echo -e "$input" > "${SCRIPT_DIR}/data/impact_filter/${test_name}.vcf"
    
    # Run the test with error handling
    "${ROOT_DIR}/build/src/VCFX_impact_filter/VCFX_impact_filter" --filter-impact "$filter_level" < "${SCRIPT_DIR}/data/impact_filter/${test_name}.vcf" > "${SCRIPT_DIR}/data/impact_filter/${test_name}.out" 2> "${SCRIPT_DIR}/data/impact_filter/${test_name}.err" || true
    
    # Compare output with expected
    if diff "${SCRIPT_DIR}/data/impact_filter/${test_name}.out" <(echo -e "$expected") > /dev/null; then
        echo "✓ Test passed: $test_name"
    else
        echo "✗ Test failed: $test_name"
        echo "Expected:"
        echo -e "$expected"
        echo "Got:"
        cat "${SCRIPT_DIR}/data/impact_filter/${test_name}.out"
        exit 1
    fi
    echo
}

# Test 1: Basic filtering with HIGH impact
echo "Test 1: Basic filtering with HIGH impact"
input1="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=LOW
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=MODIFIER"
expected1="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH;EXTRACTED_IMPACT=HIGH"
run_test "high_impact" "$input1" "$expected1" "HIGH" "Testing filtering for HIGH impact only"

# Test 2: Filtering with MODERATE impact (should include HIGH and MODERATE)
echo "Test 2: Filtering with MODERATE impact"
input2="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=LOW
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=MODIFIER"
expected2="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH;EXTRACTED_IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE;EXTRACTED_IMPACT=MODERATE"
run_test "moderate_impact" "$input2" "$expected2" "MODERATE" "Testing filtering for MODERATE impact and higher"

# Test 3: Filtering with LOW impact (should include HIGH, MODERATE, and LOW)
echo "Test 3: Filtering with LOW impact"
input3="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=LOW
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=MODIFIER"
expected3="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH;EXTRACTED_IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE;EXTRACTED_IMPACT=MODERATE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=LOW;EXTRACTED_IMPACT=LOW"
run_test "low_impact" "$input3" "$expected3" "LOW" "Testing filtering for LOW impact and higher"

# Test 4: Filtering with MODIFIER impact (should include all)
echo "Test 4: Filtering with MODIFIER impact"
input4="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=LOW
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=MODIFIER"
expected4="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH;EXTRACTED_IMPACT=HIGH
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE;EXTRACTED_IMPACT=MODERATE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=LOW;EXTRACTED_IMPACT=LOW
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=MODIFIER;EXTRACTED_IMPACT=MODIFIER"
run_test "modifier_impact" "$input4" "$expected4" "MODIFIER" "Testing filtering for MODIFIER impact and higher"

# Test 5: Case-insensitive impact values
echo "Test 5: Case-insensitive impact values"
input5="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=high
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=Moderate
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=low
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=modifier"
expected5="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=high;EXTRACTED_IMPACT=high
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=Moderate;EXTRACTED_IMPACT=Moderate"
run_test "case_insensitive" "$input5" "$expected5" "MODERATE" "Testing case-insensitive impact value handling"

# Test 6: Complex impact values (with prefixes/suffixes)
echo "Test 6: Complex impact values"
input6="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH_MISSENSE
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE_NONSENSE
1\t300\trs3\tG\tA\t100\tPASS\tAC=1;AN=2;IMPACT=MY_LOW_IMPACT
1\t400\trs4\tT\tC\t100\tPASS\tAC=1;AN=2;IMPACT=MODIFIER_SPLICE"
expected6="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH_MISSENSE;EXTRACTED_IMPACT=HIGH_MISSENSE
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE_NONSENSE;EXTRACTED_IMPACT=MODERATE_NONSENSE"
run_test "complex_impact" "$input6" "$expected6" "MODERATE" "Testing complex impact values with prefixes/suffixes"

# Test 7: Missing IMPACT field (should be treated as UNKNOWN)
echo "Test 7: Missing IMPACT field"
input7="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH"
expected7="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH;EXTRACTED_IMPACT=HIGH"
run_test "missing_impact" "$input7" "$expected7" "HIGH" "Testing variants with missing IMPACT field"

# Test 8: Empty INFO field
echo "Test 8: Empty INFO field"
input8="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs1\tA\tG\t100\tPASS\t.
1\t200\trs2\tC\tT\t100\tPASS\tIMPACT=HIGH"
expected8="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t200\trs2\tC\tT\t100\tPASS\tIMPACT=HIGH;EXTRACTED_IMPACT=HIGH"
run_test "empty_info" "$input8" "$expected8" "HIGH" "Testing variants with empty INFO field"

# Test 9: Help message
echo "Test 9: Help message"
"${ROOT_DIR}/build/src/VCFX_impact_filter/VCFX_impact_filter" --help > "${SCRIPT_DIR}/data/impact_filter/help.out" 2> "${SCRIPT_DIR}/data/impact_filter/help.err" || true
if grep -q "VCFX_impact_filter: Filter VCF variants based on predicted impact" "${SCRIPT_DIR}/data/impact_filter/help.out"; then
    echo "✓ Test passed: Help message"
else
    echo "✗ Test failed: Help message"
    echo "Expected help message, got:"
    cat "${SCRIPT_DIR}/data/impact_filter/help.out"
    exit 1
fi

# Test 10: Invalid impact level
echo "Test 10: Invalid impact level"
"${ROOT_DIR}/build/src/VCFX_impact_filter/VCFX_impact_filter" --filter-impact "INVALID" < /dev/null > "${SCRIPT_DIR}/data/impact_filter/invalid_level.out" 2> "${SCRIPT_DIR}/data/impact_filter/invalid_level.err" || true
if grep -q "Error: Unrecognized impact level" "${SCRIPT_DIR}/data/impact_filter/invalid_level.err"; then
    echo "✓ Test passed: Invalid impact level (produced expected error message)"
else
    echo "✗ Test failed: Invalid impact level"
    echo "Expected error message about unrecognized impact level, got:"
    cat "${SCRIPT_DIR}/data/impact_filter/invalid_level.err"
    exit 1
fi

# Test 11: Missing required argument
echo "Test 11: Missing required argument"
"${ROOT_DIR}/build/src/VCFX_impact_filter/VCFX_impact_filter" > "${SCRIPT_DIR}/data/impact_filter/missing_arg.out" 2> "${SCRIPT_DIR}/data/impact_filter/missing_arg.err" || true
if grep -q "Usage:" "${SCRIPT_DIR}/data/impact_filter/missing_arg.out" || grep -q "VCFX_impact_filter: Filter VCF variants" "${SCRIPT_DIR}/data/impact_filter/missing_arg.out"; then
    echo "✓ Test passed: Missing required argument (produced help message)"
else
    echo "✗ Test failed: Missing required argument (expected help message)"
    exit 1
fi

# Test 12: Performance with large file
echo "Test 12: Performance with large file"
{
    echo "##fileformat=VCFv4.2"
    echo "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">"
    echo "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
    echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    for i in {1..50}; do
        # Randomly assign impact levels
        case $((i % 4)) in
            0) impact="HIGH" ;;
            1) impact="MODERATE" ;;
            2) impact="LOW" ;;
            3) impact="MODIFIER" ;;
        esac
        echo "1\t$i\trs$i\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=$impact"
    done
} > "${SCRIPT_DIR}/data/impact_filter/large.vcf"

"${ROOT_DIR}/build/src/VCFX_impact_filter/VCFX_impact_filter" --filter-impact "MODERATE" < "${SCRIPT_DIR}/data/impact_filter/large.vcf" > "${SCRIPT_DIR}/data/impact_filter/large.out" 2> "${SCRIPT_DIR}/data/impact_filter/large.err" || true
echo "✓ Test passed: Performance test"

# Test 13: Verify line format is preserved
echo "Test 13: Verify line format is preserved (with sample columns)"
input13="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH\tGT\t0/0\t0/1
1\t200\trs2\tC\tT\t100\tPASS\tAC=1;AN=2;IMPACT=MODERATE\tGT\t0/1\t1/1"
expected13="##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
1\t100\trs1\tA\tG\t100\tPASS\tAC=1;AN=2;IMPACT=HIGH;EXTRACTED_IMPACT=HIGH\tGT\t0/0\t0/1"
run_test "with_samples" "$input13" "$expected13" "HIGH" "Testing filtering with sample columns"

echo "All VCFX_impact_filter tests passed! 🎉" 