#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_phred_filter ==="

# Paths
VCFX_PHRED_FILTER="../build/src/VCFX_phred_filter/VCFX_phred_filter"
TEST_DATA_DIR="data/phred_filter"
EXPECTED_DIR="expected/phred_filter"
OUTPUT_DIR="out/phred_filter"

# Check if executable exists
if [ ! -f "$VCFX_PHRED_FILTER" ]; then
  echo "Error: $VCFX_PHRED_FILTER not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Test counter
TOTAL_TESTS=0
PASSED_TESTS=0

# Function to run a test and compare outputs
run_test() {
  local name="$1"
  local command="$2"
  local input="$3"
  local expected="$4"
  local output="$OUTPUT_DIR/$name.out"
  
  TOTAL_TESTS=$((TOTAL_TESTS + 1))
  
  echo "Running test: $name"
  eval "$command < $input > $output"
  
  if diff -u "$expected" "$output" > /dev/null; then
    echo "✅ Test passed: $name"
    PASSED_TESTS=$((PASSED_TESTS + 1))
  else
    echo "❌ Test failed: $name"
    echo "Expected output (from $expected):"
    cat "$expected"
    echo
    echo "Actual output (from $output):"
    cat "$output"
    echo
    echo "Diff:"
    diff -u "$expected" "$output" || true
  fi
  echo
}

# Test 1: Basic filtering with threshold 30
run_test "basic_threshold_30" \
  "$VCFX_PHRED_FILTER -p 30" \
  "$TEST_DATA_DIR/basic.vcf" \
  "$EXPECTED_DIR/basic_threshold_30.vcf"

# Test 2: Basic filtering with threshold 20
run_test "basic_threshold_20" \
  "$VCFX_PHRED_FILTER -p 20" \
  "$TEST_DATA_DIR/basic.vcf" \
  "$EXPECTED_DIR/basic_threshold_20.vcf"

# Test 3: Basic filtering with threshold 30 and keep missing QUAL
run_test "basic_threshold_30_keep_missing" \
  "$VCFX_PHRED_FILTER -p 30 -k" \
  "$TEST_DATA_DIR/basic.vcf" \
  "$EXPECTED_DIR/basic_threshold_30_keep_missing.vcf"

# Test 4: Malformed QUAL values with threshold 30
run_test "malformed_threshold_30" \
  "$VCFX_PHRED_FILTER -p 30" \
  "$TEST_DATA_DIR/malformed.vcf" \
  "$EXPECTED_DIR/malformed_threshold_30.vcf"

# Test 5: Malformed QUAL values with threshold 5
run_test "malformed_threshold_5" \
  "$VCFX_PHRED_FILTER -p 5" \
  "$TEST_DATA_DIR/malformed.vcf" \
  "$EXPECTED_DIR/malformed_threshold_5.vcf"

# Test 6: Malformed QUAL values with threshold 30 and keep missing QUAL
run_test "malformed_threshold_30_keep_missing" \
  "$VCFX_PHRED_FILTER -p 30 -k" \
  "$TEST_DATA_DIR/malformed.vcf" \
  "$EXPECTED_DIR/malformed_threshold_30_keep_missing.vcf"

# Test 7: Invalid/missing records with threshold 30
run_test "invalid_records_threshold_30" \
  "$VCFX_PHRED_FILTER -p 30" \
  "$TEST_DATA_DIR/invalid_records.vcf" \
  "$EXPECTED_DIR/invalid_records_threshold_30.vcf"

# Test 8: Help message
run_test "help_message" \
  "$VCFX_PHRED_FILTER --help" \
  "/dev/null" \
  "$EXPECTED_DIR/help_message.txt"

# Test 9: Long option format
run_test "long_option_format" \
  "$VCFX_PHRED_FILTER --phred-filter 30 --keep-missing-qual" \
  "$TEST_DATA_DIR/basic.vcf" \
  "$EXPECTED_DIR/basic_threshold_30_keep_missing.vcf"

# Summary
echo "Tests completed: $PASSED_TESTS passed, $((TOTAL_TESTS - PASSED_TESTS)) failed (out of $TOTAL_TESTS total)"

if [ $PASSED_TESTS -eq $TOTAL_TESTS ]; then
  echo "✅ All VCFX_phred_filter tests passed!"
  exit 0
else
  echo "❌ Some VCFX_phred_filter tests failed."
  exit 1
fi 