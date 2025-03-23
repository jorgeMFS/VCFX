#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_variant_counter ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_variant_counter/VCFX_variant_counter"

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
      exit 1
    else
      echo "  Test $test_num passed (expected failure)"
    fi
  else
    eval "$test_cmd" > "$output_file" 2> /dev/null
    diff -u "$expected_file" "$output_file" || {
      echo "  Test $test_num failed. Actual output:"
      cat "$output_file"
      exit 1
    }
    echo "  Test $test_num passed."
  fi
}

# Create test data if it doesn't exist
if [ ! -f data/variant_counter_normal.vcf ]; then
  cat > data/variant_counter_normal.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	80	PASS	DP=60;AF=0.10	GT:DP	0/0:28	0/1:32
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05	GT:DP	0/1:15	0/0:14
chr1	40000	.	T	C	90	PASS	DP=87;AF=0.35	GT:DP	0/1:45	1/1:42
chr1	50000	rs102	G	T	30	LowQual	DP=18;AF=0.02	GT:DP	0/0:10	0/1:8
EOF
fi

if [ ! -f data/variant_counter_invalid.vcf ]; then
  cat > data/variant_counter_invalid.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
chr1	10000	rs123	A	G	100	PASS
chr1	20000	rs456	C	T	80	PASS
chr1	30000	rs789	G	A	50	LowQual
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	40000	.	T	C	90	PASS	DP=87;AF=0.35	GT:DP	0/1:45	1/1:42
chr1	50000	rs102	G	T	30	LowQual	DP=18;AF=0.02	GT:DP	0/0:10	0/1:8
EOF
fi

if [ ! -f data/variant_counter_malformed.vcf ]; then
  cat > data/variant_counter_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	invalid	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	not_a_number	PASS	DP=60;AF=0.10	GT:DP	0/0:28	0/1:32
missing_chrom_field_but_has_8_columns	30000	rs789	G	A	50	LowQual	DP=29
chr1	40000	.	T	C	90	PASS	DP=87;AF=0.35	GT:DP	0/1:45	1/1:42
EOF
fi

if [ ! -f data/variant_counter_large.vcf ]; then
  # Create a larger VCF file with 100 entries to test performance
  cat > data/variant_counter_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF
  
  # Add 100 variant lines
  for i in $(seq 1 100); do
    echo "chr1	$((10000 + $i * 100))	rs$i	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30	0/0:25" >> data/variant_counter_large.vcf
  done
fi

if [ ! -f expected/variant_counter_normal.txt ]; then
  cat > expected/variant_counter_normal.txt << EOF
Total Variants: 5
EOF
fi

if [ ! -f expected/variant_counter_invalid_nonstrict.txt ]; then
  cat > expected/variant_counter_invalid_nonstrict.txt << EOF
Total Variants: 2
EOF
fi

if [ ! -f expected/variant_counter_normal_nonstrict.txt ]; then
  cat > expected/variant_counter_normal_nonstrict.txt << EOF
Total Variants: 5
EOF
fi

if [ ! -f expected/variant_counter_empty.txt ]; then
  cat > expected/variant_counter_empty.txt << EOF
Total Variants: 0
EOF
fi

if [ ! -f expected/variant_counter_malformed_nonstrict.txt ]; then
  cat > expected/variant_counter_malformed_nonstrict.txt << EOF
Total Variants: 4
EOF
fi

if [ ! -f expected/variant_counter_large.txt ]; then
  cat > expected/variant_counter_large.txt << EOF
Total Variants: 100
EOF
fi

if [ ! -f data/variant_counter_empty.vcf ]; then
  cat > data/variant_counter_empty.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF
fi

# Test 1: Count variants in a normal VCF file (strict mode)
run_test 1 "Counting variants in a normal VCF file (strict mode)" \
  "cat data/variant_counter_normal.vcf | $VCFX_EXECUTABLE --strict" \
  "expected/variant_counter_normal.txt" \
  "out/variant_counter_normal.txt"

# Test 2: Count variants in an invalid VCF file (non-strict mode)
run_test 2 "Counting variants in an invalid VCF file (non-strict mode)" \
  "cat data/variant_counter_invalid.vcf | $VCFX_EXECUTABLE" \
  "expected/variant_counter_invalid_nonstrict.txt" \
  "out/variant_counter_invalid_nonstrict.txt"

# Test 3: Count variants in an invalid VCF file (strict mode) - should fail
run_test 3 "Counting variants in an invalid VCF file (strict mode)" \
  "cat data/variant_counter_invalid.vcf | $VCFX_EXECUTABLE --strict" \
  "/dev/null" \
  "out/variant_counter_invalid_strict.txt" \
  "true"

# Test 4: Count variants in a normal VCF file (non-strict mode)
run_test 4 "Counting variants in a normal VCF file (non-strict mode)" \
  "cat data/variant_counter_normal.vcf | $VCFX_EXECUTABLE" \
  "expected/variant_counter_normal_nonstrict.txt" \
  "out/variant_counter_normal_nonstrict.txt"

# Test 5: Count variants in an empty VCF file
run_test 5 "Counting variants in an empty VCF file" \
  "cat data/variant_counter_empty.vcf | $VCFX_EXECUTABLE" \
  "expected/variant_counter_empty.txt" \
  "out/variant_counter_empty.txt"

# Test 6: Help display
echo "Test 6: Verifying help display"
$VCFX_EXECUTABLE --help > out/variant_counter_help.txt
if ! grep -q "VCFX_variant_counter: Counts the total number of valid variants in a VCF" out/variant_counter_help.txt; then
  echo "  Test 6 failed. Help message not displayed correctly."
  cat out/variant_counter_help.txt
  exit 1
fi
echo "  Test 6 passed."

# Test 7: Malformed VCF file (non-strict mode)
run_test 7 "Counting variants in a malformed VCF file (non-strict mode)" \
  "cat data/variant_counter_malformed.vcf | $VCFX_EXECUTABLE" \
  "expected/variant_counter_malformed_nonstrict.txt" \
  "out/variant_counter_malformed_nonstrict.txt"

# Test 8: Performance with large VCF file
echo "Test 8: Performance with large VCF file"
time cat data/variant_counter_large.vcf | $VCFX_EXECUTABLE > out/variant_counter_large.txt
diff -u expected/variant_counter_large.txt out/variant_counter_large.txt || {
  echo "  Test 8 failed. Actual output:"
  cat out/variant_counter_large.txt
  exit 1
}
echo "  Test 8 passed."

echo "All VCFX_variant_counter tests passed!" 