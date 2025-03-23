#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_missing_detector ==="

# Paths
VCFX_MISSING_DETECTOR="../build/src/VCFX_missing_detector/VCFX_missing_detector"
TEST_DATA_DIR="data/missing_detector"
EXPECTED_DIR="expected/missing_detector"
OUTPUT_DIR="out/missing_detector"

# Check if executable exists
if [ ! -f "$VCFX_MISSING_DETECTOR" ]; then
  echo "Error: $VCFX_MISSING_DETECTOR not found!"
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
  
  # Remove trailing whitespace for comparison and use diff with -w -B options
  # to ignore all whitespace differences
  sed 's/[[:space:]]*$//' "$expected" > "$expected.clean"
  sed 's/[[:space:]]*$//' "$output" > "$output.clean"
  
  if diff -w -B "$expected.clean" "$output.clean" > /dev/null; then
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
    diff -w -B "$expected.clean" "$output.clean" || true
  fi
  
  # Clean up temporary files
  rm -f "$expected.clean" "$output.clean"
  echo
}

# Test 1: Basic detection of missing genotypes
run_test "basic" \
  "$VCFX_MISSING_DETECTOR" \
  "$TEST_DATA_DIR/basic.vcf" \
  "$EXPECTED_DIR/basic_out.vcf"

# Test 2: Malformed VCF handling
run_test "malformed" \
  "$VCFX_MISSING_DETECTOR" \
  "$TEST_DATA_DIR/malformed.vcf" \
  "$EXPECTED_DIR/malformed_out.vcf"

# Test 3: Empty VCF (headers only)
run_test "empty" \
  "$VCFX_MISSING_DETECTOR" \
  "$TEST_DATA_DIR/empty.vcf" \
  "$EXPECTED_DIR/empty_out.vcf"

# Test 4: Help message
run_test "help_message" \
  "$VCFX_MISSING_DETECTOR --help" \
  "/dev/null" \
  "$EXPECTED_DIR/help_message.txt"

# Summary
echo "Tests completed: $PASSED_TESTS passed, $((TOTAL_TESTS - PASSED_TESTS)) failed (out of $TOTAL_TESTS total)"

if [ $PASSED_TESTS -eq $TOTAL_TESTS ]; then
  echo "✅ All VCFX_missing_detector tests passed!"
  exit 0
else
  echo "❌ Some VCFX_missing_detector tests failed."
  exit 1
fi 