#!/usr/bin/env bash

# Track failures
failures=0

echo "=== Testing VCFX_sv_handler ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Define paths
SV_HANDLER="../build/src/VCFX_sv_handler/VCFX_sv_handler"

# Check if executable exists
if [ ! -f "$SV_HANDLER" ]; then
  echo "Error: $SV_HANDLER not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Function to run a test and report results
run_test() {
  local test_num="$1"
  local test_desc="$2"
  local test_cmd="$3"
  local expected_file="$4"
  local output_file="$5"
  local expect_failure="${6:-false}"
  
  echo "Test $test_num: $test_desc"
  
  if [ "$expect_failure" = "true" ]; then
    if eval "$test_cmd" > "$output_file" 2>&1; then
      echo "  Error: Expected failure but got success"
      ((failures++))
    else
      echo "  Test $test_num passed (expected failure)"
    fi
  else
    eval "$test_cmd" > "$output_file" 2> /dev/null
    if ! diff -u "$expected_file" "$output_file"; then
      echo "  Test $test_num failed. Actual output:"
      cat "$output_file"
      ((failures++))
    else
      echo "  Test $test_num passed."
    fi
  fi
}

# Create test data: Various SVs
if [ ! -f data/sv_test.vcf ]; then
  cat > data/sv_test.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=15000
chr1	20000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=25000
chr1	30000	.	C	<INV>	.	PASS	SVTYPE=INV;END=35000
chr1	40000	.	T	<BND>	.	PASS	SVTYPE=BND;MATEID=bnd_1
chr1	50000	.	A	G	.	PASS	.
chr1	60000	.	C	G	.	PASS	AF=0.1;DP=100
chr1	70000	.	T	<CNV>	.	PASS	SVTYPE=CNV;END=75000
chr1	80000	.	A	<DEL>	.	PASS	SVTYPE=DEL
chr1	90000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=invalid
EOF
fi

# Create expected output for filter-only
if [ ! -f expected/sv_filter_only.vcf ]; then
  cat > expected/sv_filter_only.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=15000
chr1	20000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=25000
chr1	30000	.	C	<INV>	.	PASS	SVTYPE=INV;END=35000
chr1	40000	.	T	<BND>	.	PASS	SVTYPE=BND;MATEID=bnd_1
chr1	70000	.	T	<CNV>	.	PASS	SVTYPE=CNV;END=75000
chr1	80000	.	A	<DEL>	.	PASS	SVTYPE=DEL
chr1	90000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=invalid
EOF
fi

# Create expected output for sv-modify
if [ ! -f expected/sv_modify.vcf ]; then
  cat > expected/sv_modify.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=15000;SV_VALIDATED=1;SV_SIZE=5000
chr1	20000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=25000;SV_VALIDATED=1;SV_SIZE=5000
chr1	30000	.	C	<INV>	.	PASS	SVTYPE=INV;END=35000;SV_VALIDATED=1;INV_TYPE=PARALLEL
chr1	40000	.	T	<BND>	.	PASS	SVTYPE=BND;MATEID=bnd_1;SV_VALIDATED=1;BND_ORIENTATION=PAIR
chr1	50000	.	A	G	.	PASS	.
chr1	60000	.	C	G	.	PASS	AF=0.1;DP=100
chr1	70000	.	T	<CNV>	.	PASS	SVTYPE=CNV;END=75000;SV_VALIDATED=1
chr1	80000	.	A	<DEL>	.	PASS	SVTYPE=DEL;SV_VALIDATED=1
chr1	90000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=invalid;SV_VALIDATED=1
EOF
fi

# Create expected output for both filter and modify
if [ ! -f expected/sv_filter_modify.vcf ]; then
  cat > expected/sv_filter_modify.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=15000;SV_VALIDATED=1;SV_SIZE=5000
chr1	20000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=25000;SV_VALIDATED=1;SV_SIZE=5000
chr1	30000	.	C	<INV>	.	PASS	SVTYPE=INV;END=35000;SV_VALIDATED=1;INV_TYPE=PARALLEL
chr1	40000	.	T	<BND>	.	PASS	SVTYPE=BND;MATEID=bnd_1;SV_VALIDATED=1;BND_ORIENTATION=PAIR
chr1	70000	.	T	<CNV>	.	PASS	SVTYPE=CNV;END=75000;SV_VALIDATED=1
chr1	80000	.	A	<DEL>	.	PASS	SVTYPE=DEL;SV_VALIDATED=1
chr1	90000	.	G	<DUP>	.	PASS	SVTYPE=DUP;END=invalid;SV_VALIDATED=1
EOF
fi

# Create test data with invalid SVs
if [ ! -f data/sv_invalid.vcf ]; then
  cat > data/sv_invalid.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	invalid	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=15000
chr1	20000	.	G	<DUP>	.	PASS	SVTYPE=
chr1	30000	.	AAAA	<INV>	.	PASS	SVTYPE=INV;END=35000
EOF
fi

# Run the tests
# Test 1: Filter only (--sv-filter-only)
run_test 1 "Filter only structural variants" \
  "$SV_HANDLER --sv-filter-only < data/sv_test.vcf" \
  "expected/sv_filter_only.vcf" \
  "out/sv_filter_only.vcf"

# Test 2: Modify SVs (--sv-modify)
run_test 2 "Modify structural variants" \
  "$SV_HANDLER --sv-modify < data/sv_test.vcf" \
  "expected/sv_modify.vcf" \
  "out/sv_modify.vcf"

# Test 3: Both filter and modify
run_test 3 "Filter and modify structural variants" \
  "$SV_HANDLER --sv-filter-only --sv-modify < data/sv_test.vcf" \
  "expected/sv_filter_modify.vcf" \
  "out/sv_filter_modify.vcf"

# Test 4: Help message
echo "Test 4: Help message"
$SV_HANDLER --help > out/sv_help.txt
if ! grep -q "VCFX_sv_handler" out/sv_help.txt || ! grep -q "sv-filter-only" out/sv_help.txt; then
  echo "  Test 4 failed: Help message doesn't contain expected content"
  ((failures++))
else
  echo "  Test 4 passed."
fi

# Test 5: Invalid inputs (expect warnings but not crashes)
echo "Test 5: Handling invalid inputs"
$SV_HANDLER --sv-modify < data/sv_invalid.vcf > out/sv_invalid_out.vcf 2> out/sv_invalid_warnings.txt
if [ $? -ne 0 ]; then
  echo "  Test 5 failed: Program crashed on invalid input"
  ((failures++))
else
  if ! grep -q "Warning" out/sv_invalid_warnings.txt; then
    echo "  Test 5 failed: Expected warnings for invalid input"
    ((failures++))
  else
    echo "  Test 5 passed: Program handled invalid input gracefully"
  fi
fi

# Report results
if [ $failures -eq 0 ]; then
  echo "All VCFX_sv_handler tests passed!"
  exit 0
else
  echo "$failures VCFX_sv_handler tests failed."
  exit 1
fi 