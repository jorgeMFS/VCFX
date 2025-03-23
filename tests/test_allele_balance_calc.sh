#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_allele_balance_calc ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_allele_balance_calc/VCFX_allele_balance_calc"

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
$VCFX_EXECUTABLE \
  < data/allele_balance_calc_A.vcf \
  > out/allele_balance_calc_A_out.tsv

diff -u expected/allele_balance_calc_A_out.tsv out/allele_balance_calc_A_out.tsv
echo "Scenario A passed."

# Scenario B
$VCFX_EXECUTABLE \
  < data/allele_balance_calc_B.vcf \
  > out/allele_balance_calc_B_out.tsv

diff -u expected/allele_balance_calc_B_out.tsv out/allele_balance_calc_B_out.tsv
echo "Scenario B passed."

echo "All allele_balance_calc tests passed!"
