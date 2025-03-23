#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_allele_balance_filter ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_allele_balance_filter/VCFX_allele_balance_filter"

# Check if executable exists
if [ ! -f "$VCFX_EXECUTABLE" ]; then
  echo "Error: $VCFX_EXECUTABLE not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Create output directories if they don't exist
mkdir -p out
mkdir -p data
mkdir -p expected

# Scenario A
$VCFX_EXECUTABLE --filter-allele-balance 0.5 \
  < data/allele_balance_filter_A.vcf \
  > out/allele_balance_filter_A_out.vcf

diff -u expected/allele_balance_filter_A_out.vcf out/allele_balance_filter_A_out.vcf
echo "Scenario A passed."

# Scenario B
$VCFX_EXECUTABLE --filter-allele-balance 0.1 \
  < data/allele_balance_filter_B.vcf \
  > out/allele_balance_filter_B_out.vcf

diff -u expected/allele_balance_filter_B_out.vcf out/allele_balance_filter_B_out.vcf
echo "Scenario B passed."

echo "All allele_balance_filter tests passed!"
