#!/usr/bin/env bash

# VCFX Population Filter Tests
# Test that the VCFX_population_filter tool correctly filters VCF files by population

# Exit on error
set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
BUILD_DIR="../build"
TOOL="${BUILD_DIR}/src/VCFX_population_filter/VCFX_population_filter"  # Fixed path based on actual location
DATA_DIR="data/population_filter"
EXPECTED_DIR="expected/population_filter"
TMP_DIR="tmp"

# Check if build exists and compile if not
if [ ! -e "$TOOL" ]; then
    echo "VCFX_population_filter not found. Building first..."
    cd ..
    mkdir -p build
    cd build
    cmake ..
    make VCFX_population_filter -j
    cd "$SCRIPT_DIR"
fi

# Prepare tmp directory
mkdir -p "$TMP_DIR"

# Test 1: Filter sample.vcf to include only EUR samples
echo "Test 1: EUR Population Filter"
"$TOOL" --population EUR --pop-map "$DATA_DIR/pop_map.txt" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/eur_output.vcf"
if diff -q "$TMP_DIR/eur_output.vcf" "$EXPECTED_DIR/eur_filter.vcf" > /dev/null; then
    echo "✅ Test 1 passed: EUR population filtering works correctly"
else
    echo "❌ Test 1 failed: Unexpected output for EUR population filtering"
    echo "Expected:"
    cat "$EXPECTED_DIR/eur_filter.vcf"
    echo "Got:"
    cat "$TMP_DIR/eur_output.vcf"
    exit 1
fi

# Test 2: Filter sample.vcf to include only AFR samples
echo "Test 2: AFR Population Filter"
"$TOOL" --population AFR --pop-map "$DATA_DIR/pop_map.txt" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/afr_output.vcf"
if diff -q "$TMP_DIR/afr_output.vcf" "$EXPECTED_DIR/afr_filter.vcf" > /dev/null; then
    echo "✅ Test 2 passed: AFR population filtering works correctly"
else
    echo "❌ Test 2 failed: Unexpected output for AFR population filtering"
    echo "Expected:"
    cat "$EXPECTED_DIR/afr_filter.vcf"
    echo "Got:"
    cat "$TMP_DIR/afr_output.vcf"
    exit 1
fi

# Test 3: Filter sample.vcf to include only EAS samples
echo "Test 3: EAS Population Filter"
"$TOOL" --population EAS --pop-map "$DATA_DIR/pop_map.txt" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/eas_output.vcf"
if diff -q "$TMP_DIR/eas_output.vcf" "$EXPECTED_DIR/eas_filter.vcf" > /dev/null; then
    echo "✅ Test 3 passed: EAS population filtering works correctly"
else
    echo "❌ Test 3 failed: Unexpected output for EAS population filtering"
    echo "Expected:"
    cat "$EXPECTED_DIR/eas_filter.vcf"
    echo "Got:"
    cat "$TMP_DIR/eas_output.vcf"
    exit 1
fi

# Test 4: Filter sample.vcf with unknown population
echo "Test 4: Unknown Population Filter"
"$TOOL" --population XYZ --pop-map "$DATA_DIR/pop_map.txt" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/unknown_output.vcf"
if diff -q "$TMP_DIR/unknown_output.vcf" "$EXPECTED_DIR/unknown_filter.vcf" > /dev/null; then
    echo "✅ Test 4 passed: Unknown population filtering works correctly"
else
    echo "❌ Test 4 failed: Unexpected output for unknown population filtering"
    echo "Expected:"
    cat "$EXPECTED_DIR/unknown_filter.vcf"
    echo "Got:"
    cat "$TMP_DIR/unknown_output.vcf"
    exit 1
fi

# Test 5: Filter malformed.vcf file
echo "Test 5: Malformed VCF handling"
"$TOOL" --population EUR --pop-map "$DATA_DIR/pop_map.txt" < "$DATA_DIR/malformed.vcf" > "$TMP_DIR/malformed_output.vcf"
if diff -q "$TMP_DIR/malformed_output.vcf" "$EXPECTED_DIR/malformed_filter.vcf" > /dev/null; then
    echo "✅ Test 5 passed: Malformed VCF handling works correctly"
else
    echo "❌ Test 5 failed: Unexpected output for malformed VCF handling"
    echo "Expected:"
    cat "$EXPECTED_DIR/malformed_filter.vcf"
    echo "Got:"
    cat "$TMP_DIR/malformed_output.vcf"
    exit 1
fi

# Test 6: Check help message
echo "Test 6: Help message"
"$TOOL" --help > "$TMP_DIR/help_output.txt"
if grep -q "VCFX_population_filter: Subset VCF to samples in specified population" "$TMP_DIR/help_output.txt"; then
    echo "✅ Test 6 passed: Help message displayed correctly"
else
    echo "❌ Test 6 failed: Help message not displayed correctly"
    echo "Got:"
    cat "$TMP_DIR/help_output.txt"
    exit 1
fi

# Test 7: No arguments should display help message
echo "Test 7: No arguments handling"
"$TOOL" > "$TMP_DIR/no_args_output.txt"
if grep -q "VCFX_population_filter: Subset VCF to samples in specified population" "$TMP_DIR/no_args_output.txt"; then
    echo "✅ Test 7 passed: Help message displayed when no arguments provided"
else
    echo "❌ Test 7 failed: Help message not displayed when no arguments provided"
    echo "Got:"
    cat "$TMP_DIR/no_args_output.txt"
    exit 1
fi

# Test 8: Empty population map handling
echo "Test 8: Empty population map handling"
if "$TOOL" --population EUR --pop-map "$DATA_DIR/empty_pop_map.txt" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/empty_map_output.vcf" 2> "$TMP_DIR/empty_map_error.txt"; then
    if grep -q "Warning: No samples found for population tag: EUR" "$TMP_DIR/empty_map_error.txt"; then
        echo "✅ Test 8 passed: Empty population map warning displayed correctly"
    else
        echo "❌ Test 8 failed: Empty population map warning not displayed correctly"
        echo "Got stderr:"
        cat "$TMP_DIR/empty_map_error.txt"
        exit 1
    fi
else
    echo "❌ Test 8 failed: Tool should not exit with error on empty population map"
    exit 1
fi

# If we get here, all tests passed
echo "All population filter tests passed! 🎉"
exit 0 