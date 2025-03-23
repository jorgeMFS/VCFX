#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_distance_calculator ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_distance_calculator/VCFX_distance_calculator"

# Check if executable exists
if [ ! -f "$VCFX_EXECUTABLE" ]; then
  echo "Error: $VCFX_EXECUTABLE not found!"
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
if [ ! -f data/distance_normal.vcf ]; then
  cat > data/distance_normal.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	80	PASS	DP=60;AF=0.10	GT:DP	0/0:28	0/1:32
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05	GT:DP	0/1:15	0/0:14
chr2	5000	rs111	T	C	90	PASS	DP=87;AF=0.35	GT:DP	0/1:45	1/1:42
chr2	15000	rs222	G	T	30	LowQual	DP=18;AF=0.02	GT:DP	0/0:10	0/1:8
chr2	25000	rs333	A	C	60	PASS	DP=40;AF=0.15	GT:DP	0/1:20	0/1:20
EOF
fi

# Create multi-chromosome test data
if [ ! -f data/distance_multi_chrom.vcf ]; then
  cat > data/distance_multi_chrom.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	rs1	A	G	100	PASS	.
chr2	15000	rs2	C	T	90	PASS	.
chr1	25000	rs3	G	A	80	PASS	.
chr3	5000	rs4	T	C	70	PASS	.
chr2	35000	rs5	G	T	60	PASS	.
chr3	15000	rs6	A	C	50	PASS	.
chr1	45000	rs7	C	G	40	PASS	.
EOF
fi

# Create sorted test data (for testing consecutively-ordered variants)
if [ ! -f data/distance_sorted.vcf ]; then
  cat > data/distance_sorted.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	rs1	A	G	100	PASS	.
chr1	2000	rs2	C	T	90	PASS	.
chr1	4000	rs3	G	A	80	PASS	.
chr1	8000	rs4	T	C	70	PASS	.
chr1	16000	rs5	G	T	60	PASS	.
EOF
fi

# Create unsorted test data (for testing that variants don't need to be pre-sorted)
if [ ! -f data/distance_unsorted.vcf ]; then
  cat > data/distance_unsorted.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	16000	rs5	G	T	60	PASS	.
chr1	4000	rs3	G	A	80	PASS	.
chr1	8000	rs4	T	C	70	PASS	.
chr1	1000	rs1	A	G	100	PASS	.
chr1	2000	rs2	C	T	90	PASS	.
EOF
fi

# Create adjacent variants test data (distance = 1)
if [ ! -f data/distance_adjacent.vcf ]; then
  cat > data/distance_adjacent.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	rs1	A	G	100	PASS	.
chr1	1001	rs2	C	T	90	PASS	.
chr1	1002	rs3	G	A	80	PASS	.
chr1	1003	rs4	T	C	70	PASS	.
chr2	2000	rs5	G	T	60	PASS	.
chr2	2001	rs6	A	C	50	PASS	.
EOF
fi

# Create test data with malformed lines
if [ ! -f data/distance_malformed.vcf ]; then
  cat > data/distance_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	rs1	A	G	100	PASS	.
not_a_chromosome	20000	rs2	C	T	90	PASS	.
chr2	not_a_position	rs3	G	A	80	PASS	.
chr2	30000	rs4	T	C	70	PASS	.
malformed_line_with_insufficient_columns
chr3	40000	rs5	G	T	60	PASS	.
EOF
fi

# Create test data with no header
if [ ! -f data/distance_no_header.vcf ]; then
  cat > data/distance_no_header.vcf << EOF
chr1	10000	rs1	A	G	100	PASS	.
chr1	20000	rs2	C	T	90	PASS	.
chr1	30000	rs3	G	A	80	PASS	.
EOF
fi

# Create test data with only one variant per chromosome
if [ ! -f data/distance_single_per_chrom.vcf ]; then
  cat > data/distance_single_per_chrom.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	rs1	A	G	100	PASS	.
chr2	20000	rs2	C	T	90	PASS	.
chr3	30000	rs3	G	A	80	PASS	.
EOF
fi

# Create test data for performance testing (many variants)
if [ ! -f data/distance_large.vcf ]; then
  cat > data/distance_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
  
  # Add 1000 variant lines with increasing positions
  for i in $(seq 1 1000); do
    echo "chr$((1 + i % 5))\t$((10000 + i * 100))\trs$i\tA\tG\t100\tPASS\t." >> data/distance_large.vcf
  done
fi

# Create test data with zero distance (duplicate positions)
if [ ! -f data/distance_duplicates.vcf ]; then
  cat > data/distance_duplicates.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	rs1	A	G	100	PASS	.
chr1	10000	rs2	C	T	90	PASS	.
chr1	20000	rs3	G	A	80	PASS	.
chr1	20000	rs4	T	C	70	PASS	.
chr2	15000	rs5	G	T	60	PASS	.
chr2	15000	rs6	A	C	50	PASS	.
EOF
fi

# Expected outputs
if [ ! -f expected/distance_normal.tsv ]; then
  cat > expected/distance_normal.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	10000	NA	NA
chr1	20000	10000	10000
chr1	30000	20000	10000
chr2	5000	NA	NA
chr2	15000	5000	10000
chr2	25000	15000	10000
EOF
fi

if [ ! -f expected/distance_multi_chrom.tsv ]; then
  cat > expected/distance_multi_chrom.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	10000	NA	NA
chr2	15000	NA	NA
chr1	25000	10000	15000
chr3	5000	NA	NA
chr2	35000	15000	20000
chr3	15000	5000	10000
chr1	45000	25000	20000
EOF
fi

if [ ! -f expected/distance_sorted.tsv ]; then
  cat > expected/distance_sorted.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	1000	NA	NA
chr1	2000	1000	1000
chr1	4000	2000	2000
chr1	8000	4000	4000
chr1	16000	8000	8000
EOF
fi

if [ ! -f expected/distance_unsorted.tsv ]; then
  cat > expected/distance_unsorted.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	16000	NA	NA
chr1	4000	16000	-12000
chr1	8000	4000	4000
chr1	1000	8000	-7000
chr1	2000	1000	1000
EOF
fi

if [ ! -f expected/distance_adjacent.tsv ]; then
  cat > expected/distance_adjacent.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	1000	NA	NA
chr1	1001	1000	1
chr1	1002	1001	1
chr1	1003	1002	1
chr2	2000	NA	NA
chr2	2001	2000	1
EOF
fi

if [ ! -f expected/distance_malformed.tsv ]; then
  cat > expected/distance_malformed.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	10000	NA	NA
chr2	30000	NA	NA
chr3	40000	NA	NA
EOF
fi

if [ ! -f expected/distance_single_per_chrom.tsv ]; then
  cat > expected/distance_single_per_chrom.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	10000	NA	NA
chr2	20000	NA	NA
chr3	30000	NA	NA
EOF
fi

if [ ! -f expected/distance_duplicates.tsv ]; then
  cat > expected/distance_duplicates.tsv << EOF
CHROM	POS	PREV_POS	DISTANCE
chr1	10000	NA	NA
chr1	10000	10000	0
chr1	20000	10000	10000
chr1	20000	20000	0
chr2	15000	NA	NA
chr2	15000	15000	0
EOF
fi

# Test 1: Normal VCF with multiple chromosomes
run_test 1 "Normal VCF with multiple chromosomes" \
  "$VCFX_EXECUTABLE < data/distance_normal.vcf" \
  "expected/distance_normal.tsv" \
  "out/distance_normal.tsv"

# Test 2: Mixed chromosome order
run_test 2 "Mixed chromosome order" \
  "$VCFX_EXECUTABLE < data/distance_multi_chrom.vcf" \
  "expected/distance_multi_chrom.tsv" \
  "out/distance_multi_chrom.tsv"

# Test 3: Sequentially sorted positions
run_test 3 "Sequentially sorted positions" \
  "$VCFX_EXECUTABLE < data/distance_sorted.vcf" \
  "expected/distance_sorted.tsv" \
  "out/distance_sorted.tsv"

# Test 4: Unsorted positions (to verify handling of non-ordered input)
run_test 4 "Unsorted positions" \
  "$VCFX_EXECUTABLE < data/distance_unsorted.vcf" \
  "expected/distance_unsorted.tsv" \
  "out/distance_unsorted.tsv"

# Test 5: Adjacent variants (distance=1)
run_test 5 "Adjacent variants (distance=1)" \
  "$VCFX_EXECUTABLE < data/distance_adjacent.vcf" \
  "expected/distance_adjacent.tsv" \
  "out/distance_adjacent.tsv"

# Test 6: Malformed input (should skip invalid lines)
run_test 6 "Malformed input (should skip invalid lines)" \
  "$VCFX_EXECUTABLE < data/distance_malformed.vcf" \
  "expected/distance_malformed.tsv" \
  "out/distance_malformed.tsv"

# Test 7: No header (should fail)
run_test 7 "No header (should fail)" \
  "$VCFX_EXECUTABLE < data/distance_no_header.vcf" \
  "/dev/null" \
  "out/distance_no_header.tsv" \
  "true"

# Test 8: Single variant per chromosome
run_test 8 "Single variant per chromosome" \
  "$VCFX_EXECUTABLE < data/distance_single_per_chrom.vcf" \
  "expected/distance_single_per_chrom.tsv" \
  "out/distance_single_per_chrom.tsv"

# Test 9: Duplicate positions (distance=0)
run_test 9 "Duplicate positions (distance=0)" \
  "$VCFX_EXECUTABLE < data/distance_duplicates.vcf" \
  "expected/distance_duplicates.tsv" \
  "out/distance_duplicates.tsv"

# Test 10: Performance with large file
echo "Test 10: Performance with large file"
time $VCFX_EXECUTABLE < data/distance_large.vcf > out/distance_large.tsv
# Just check that output is non-empty and contains expected format
if ! grep -q "CHROM" out/distance_large.tsv || ! grep -q "DISTANCE" out/distance_large.tsv; then
  echo "  Test 10 failed. Output doesn't contain expected headers."
  exit 1
fi

# Count the number of data lines (excluding header)
data_lines=$(grep -v "CHROM" out/distance_large.tsv | wc -l)
if [ "$data_lines" -lt 990 ]; then  # Should have ~1000 lines
  echo "  Test 10 failed. Expected at least 990 data lines, but found $data_lines."
  exit 1
fi
echo "  Test 10 passed."

# Test 11: Help display
echo "Test 11: Verifying help display"
$VCFX_EXECUTABLE --help > out/distance_help.txt
if ! grep -q "VCFX_distance_calculator" out/distance_help.txt; then
  echo "  Test 11 failed. Help message not displayed correctly."
  cat out/distance_help.txt
  exit 1
fi
echo "  Test 11 passed."

echo "All VCFX_distance_calculator tests passed!" 