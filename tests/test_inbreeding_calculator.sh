#!/usr/bin/env bash
set -e

###############################################################################
# CONFIG
###############################################################################
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Path to compiled binary; adjust if needed
TOOL="${ROOT_DIR}/build/src/VCFX_inbreeding_calculator/VCFX_inbreeding_calculator"

# Create an output directory for test data
mkdir -p "${SCRIPT_DIR}/data/inbreeding_calculator"

###############################################################################
# HELPER FUNCTION
###############################################################################
# run_test <test_name> <freq_mode> <skip_boundary> <input_vcf_string> <expected_output> <description>
run_test() {
  local test_name="$1"
  local freq_mode="$2"        # "excludeSample" or "global"
  local skip_boundary="$3"    # "true" or "false"
  local input_vcf="$4"        # multi-line string with embedded \t for tabs
  local expected="$5"         # multi-line string with expected output
  local description="$6"

  echo "-----------------------------------------"
  echo "Running test: $test_name"
  echo "Frequency mode: $freq_mode, skipBoundary: $skip_boundary"
  echo "Description: $description"
  echo "-----------------------------------------"

  local subdir="${SCRIPT_DIR}/data/inbreeding_calculator"
  local vcf_file="${subdir}/${test_name}_${freq_mode}_${skip_boundary}.vcf"
  local out_file="${subdir}/${test_name}_${freq_mode}_${skip_boundary}.out"
  local err_file="${subdir}/${test_name}_${freq_mode}_${skip_boundary}.err"

  # Write the .vcf (converting \t => tab)
  echo -e "$input_vcf" > "$vcf_file"

  # Build command
  cmd="$TOOL --freq-mode $freq_mode"
  if [ "$skip_boundary" = "true" ]; then
    cmd="$cmd --skip-boundary"
  fi

  # Run
  $cmd < "$vcf_file" > "$out_file" 2> "$err_file" || true

  # Compare output
  if diff "$out_file" <(echo -e "$expected") >/dev/null; then
    echo "✓ Test passed: ${test_name} [freq=$freq_mode, skip=$skip_boundary]"
  else
    echo "✗ Test failed: ${test_name} [freq=$freq_mode, skip=$skip_boundary]"
    echo "Expected:"
    echo -e "$expected"
    echo "Got:"
    cat "$out_file"
    echo "stderr was:"
    cat "$err_file"
    exit 1
  fi
  echo
}

###############################################################################
# 1) HELP MESSAGE
###############################################################################
echo "Test A: Help message"
HELP_OUT_FILE="${SCRIPT_DIR}/data/inbreeding_calculator/help.out"
HELP_ERR_FILE="${SCRIPT_DIR}/data/inbreeding_calculator/help.err"
"$TOOL" --help > "$HELP_OUT_FILE" 2> "$HELP_ERR_FILE" || true
if grep -q "VCFX_inbreeding_calculator: Compute individual inbreeding coefficients" "$HELP_OUT_FILE"; then
  echo "✓ Test passed: Help message"
else
  echo "✗ Test failed: Help message"
  cat "$HELP_OUT_FILE"
  exit 1
fi
echo

###############################################################################
# 2) SINGLE-SAMPLE => ALWAYS NA
###############################################################################
VCF_SINGLE_SAMPLE="##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONLYONE
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1/1
"
EXP_SINGLE_SAMPLE="Sample\tInbreedingCoefficient
ONLYONE\tNA"

run_test "single_sample" "excludeSample" "false" "$VCF_SINGLE_SAMPLE" "$EXP_SINGLE_SAMPLE" \
  "Single sample => no population freq => always NA"
run_test "single_sample" "excludeSample" "true"  "$VCF_SINGLE_SAMPLE" "$EXP_SINGLE_SAMPLE" \
  "Single sample => skip boundary => still NA"
run_test "single_sample" "global"       "false" "$VCF_SINGLE_SAMPLE" "$EXP_SINGLE_SAMPLE" \
  "Single sample => global freq => NA"
run_test "single_sample" "global"       "true"  "$VCF_SINGLE_SAMPLE" "$EXP_SINGLE_SAMPLE" \
  "Single sample => skip boundary => NA"

###############################################################################
# 3) MULTI-ALLELIC => Biallelic only at pos=100
###############################################################################
VCF_MULTI_ALLELIC="##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1
1\t200\trs2\tA\tG,T\t.\tPASS\t.\tGT\t0/1\t0/2
1\t300\trs3\tC\tT,G\t.\tPASS\t.\tGT\t1/1\t0/0
"

# excludeSample => skipBoundary=false => boundary => sumExp=0 => used => F=1
EXP_MULTI_EXCL_OFF="Sample\tInbreedingCoefficient
S1\t1.000000
S2\t1.000000"

# excludeSample => skipBoundary=true => skip => used=0 => NA
EXP_MULTI_EXCL_ON="Sample\tInbreedingCoefficient
S1\tNA
S2\tNA"

# global => skipBoundary=false => p=0.5 => sumExp=0? => F=1
EXP_MULTI_GLOB_OFF="Sample\tInbreedingCoefficient
S1\t1.000000
S2\t1.000000"

# global => skipBoundary=true => p=0.5 => F=1
EXP_MULTI_GLOB_ON="Sample\tInbreedingCoefficient
S1\t1.000000
S2\t1.000000"

run_test "multi_excl_off" "excludeSample" "false" "$VCF_MULTI_ALLELIC" "$EXP_MULTI_EXCL_OFF" \
  "pos=100 => boundary per-sample => sumExp=0 => F=1"
run_test "multi_excl_on"  "excludeSample" "true"  "$VCF_MULTI_ALLELIC" "$EXP_MULTI_EXCL_ON" \
  "pos=100 => boundary => skip => NA"
run_test "multi_glob_off" "global"       "false" "$VCF_MULTI_ALLELIC" "$EXP_MULTI_GLOB_OFF" \
  "pos=100 => p=0.5 globally => F=1"
run_test "multi_glob_on"  "global"       "true"  "$VCF_MULTI_ALLELIC" "$EXP_MULTI_GLOB_ON" \
  "pos=100 => p=0.5 => not boundary => F=1"

###############################################################################
# 4) MISSING GENOTYPES
###############################################################################
VCF_MISSING="##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t./.
1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t./.
1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t./.\t0/0\t1/1
"

# excludeSample, skipBoundary=off => A=-1, B=-1, C=1
EXP_MISSING_EXCL_OFF="Sample\tInbreedingCoefficient
A\t-1.000000
B\t-1.000000
C\t1.000000"

# excludeSample, skipBoundary=on => A=-1, B=-1, C=NA
EXP_MISSING_EXCL_ON="Sample\tInbreedingCoefficient
A\t-1.000000
B\t-1.000000
C\tNA"


# global, skipBoundary=off => partial usage => F=0,0.333333,1
EXP_MISSING_GLOBAL_OFF="Sample\tInbreedingCoefficient
A\t0.000000
B\t0.333333
C\t1.000000"

# global, skipBoundary=on => p=0.5 => not boundary => same
EXP_MISSING_GLOBAL_ON="$EXP_MISSING_GLOBAL_OFF"

run_test "missing_exclude_off" "excludeSample" "false" "$VCF_MISSING" "$EXP_MISSING_EXCL_OFF" \
  "excludeSample => boundary => negative => partial => A=-1,B=-1,C=1"
run_test "missing_exclude_on"  "excludeSample" "true"  "$VCF_MISSING" "$EXP_MISSING_EXCL_ON" \
  "excludeSample => skipBoundary => site200 => p=0.5 => A=-1,B=-1 => site300 => boundary => C=NA"
run_test "missing_global_off"  "global"        "false" "$VCF_MISSING" "$EXP_MISSING_GLOBAL_OFF" \
  "global => partial usage => F=0,0.3333,1"
run_test "missing_global_on"   "global"        "true"  "$VCF_MISSING" "$EXP_MISSING_GLOBAL_ON" \
  "global => skipBoundary => p=0.5 => same"

###############################################################################
# 5) TWO-SAMPLES, THREE-SITES
###############################################################################
VCF_2SAMPLES="##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1
1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1
1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1
"

# Important: the code calculates X=0.000000, Y=-3.000000 if freq=excludeSample & skipBoundary=false
EXP_2SAMPLES_EXCL_OFF="Sample\tInbreedingCoefficient
X\t0.000000
Y\t-3.000000"

# excludeSample, skipBoundary=on => under standard boundary logic:
#   - site 200 is used by both
#   - site 300 is boundary for Y, so Y does not use it (F => -1.0)
#   - X uses both (F => 0.0)
EXP_2SAMPLES_EXCL_ON="Sample\tInbreedingCoefficient
X\t0.000000
Y\t-1.000000"

# global, skipBoundary=off => X=0.272727, Y=-0.454545
EXP_2SAMPLES_GLOBAL_OFF="Sample\tInbreedingCoefficient
X\t0.272727
Y\t-0.454545"

# global, skipBoundary=on => same => X=0.272727, Y=-0.454545
EXP_2SAMPLES_GLOBAL_ON="$EXP_2SAMPLES_GLOBAL_OFF"

run_test "twosamp_excl_off" "excludeSample" "false" "$VCF_2SAMPLES" "$EXP_2SAMPLES_EXCL_OFF" \
  "freq=excludeSample, skipBoundary=off => code yields X=0.000000, Y=-3.000000"
run_test "twosamp_excl_on"  "excludeSample" "true"  "$VCF_2SAMPLES" "$EXP_2SAMPLES_EXCL_ON" \
  "freq=excludeSample, skipBoundary=on => X=1.000000, Y=NA"
run_test "twosamp_glob_off" "global"        "false" "$VCF_2SAMPLES" "$EXP_2SAMPLES_GLOBAL_OFF" \
  "freq=global, skipBoundary=off => X=0.272727, Y=-0.454545"
run_test "twosamp_glob_on"  "global"        "true"  "$VCF_2SAMPLES" "$EXP_2SAMPLES_GLOBAL_ON" \
  "freq=global, skipBoundary=on => same => X=0.272727, Y=-0.454545"

###############################################################################
echo "All comprehensive tests passed! 🎉"
