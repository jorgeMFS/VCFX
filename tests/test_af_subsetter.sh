#!/usr/bin/env bash
set -e

# We assume you run this script from tests/ or from the root with: tests/test_af_subsetter.sh

echo "=== Testing VCFX_af_subsetter ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_af_subsetter/VCFX_af_subsetter"

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

# Scenario A: Single-allelic
echo "[A] Single-allelic range 0.01-0.1"

$VCFX_EXECUTABLE --af-filter 0.01-0.1 \
  < data/af_subsetter_A.vcf \
  > out/af_subsetter_A_out.vcf

diff -u expected/af_subsetter_A_out.vcf out/af_subsetter_A_out.vcf
echo "  Scenario A passed."

# Scenario B: Multi-allelic with comma
echo "[B] Multi-allelic range 0.01-0.05"

$VCFX_EXECUTABLE --af-filter 0.01-0.05 \
  < data/af_subsetter_B.vcf \
  > out/af_subsetter_B_out.vcf

diff -u expected/af_subsetter_B_out.vcf out/af_subsetter_B_out.vcf
echo "  Scenario B passed."

# Scenario C: Missing or no AF
echo "[C] Missing or no AF => skip"

$VCFX_EXECUTABLE --af-filter 0.01-0.05 \
  < data/af_subsetter_C.vcf \
  > out/af_subsetter_C_out.vcf

diff -u expected/af_subsetter_C_out.vcf out/af_subsetter_C_out.vcf
echo "  Scenario C passed."

# Scenario D: Boundary test
echo "[D] Full range 0.0-1.0 => keep all"

$VCFX_EXECUTABLE --af-filter 0.0-1.0 \
  < data/af_subsetter_D.vcf \
  > out/af_subsetter_D_out.vcf

diff -u data/af_subsetter_D.vcf out/af_subsetter_D_out.vcf \
  || (echo "  Scenario D: boundary test failed" && exit 1)
echo "  Scenario D passed."

echo "All VCFX_af_subsetter scenarios passed!"
