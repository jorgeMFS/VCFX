#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_alignment_checker ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_alignment_checker/VCFX_alignment_checker"

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

#######################################
# Scenario X: Perfectly matching references
#######################################
echo "[X] Perfect reference match"

$VCFX_EXECUTABLE --alignment-discrepancy data/align_X.vcf data/align_refX.fa \
  > out/align_X.txt

diff -u expected/align_X.txt out/align_X.txt
echo "  Scenario X passed."

#######################################
# Scenario Y: Real mismatches
#######################################
echo "[Y] Mismatched bases, multi-allelic"

$VCFX_EXECUTABLE --alignment-discrepancy data/align_Y.vcf data/align_refY.fa \
  > out/align_Y.txt

diff -u expected/align_Y.txt out/align_Y.txt
echo "  Scenario Y passed."

#######################################
# Scenario Z: Out of range positions
#######################################
echo "[Z] Out-of-range positions => warnings"

$VCFX_EXECUTABLE --alignment-discrepancy data/align_Z.vcf data/align_refZ.fa \
  > out/align_Z.txt

diff -u expected/align_Z.txt out/align_Z.txt
echo "  Scenario Z passed."

echo "All alignment_checker scenarios passed!"
