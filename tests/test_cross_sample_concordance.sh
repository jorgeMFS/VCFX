#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_cross_sample_concordance ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_cross_sample_concordance/VCFX_cross_sample_concordance"

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
if [ ! -f data/concordance_all_match.vcf ]; then
  cat > data/concordance_all_match.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0/1:30:99	0/1:25:95	0/1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0/0:28:99	0/0:32:98	0/0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1/1:15:45	1/1:14:52	1/1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0/1:45:99	1/0:42:99	0/1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1/1:10:35	1/1:8:30	1/1:9:32
EOF
fi

if [ ! -f data/concordance_some_mismatch.vcf ]; then
  cat > data/concordance_some_mismatch.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0/1:30:99	0/1:25:95	0/1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0/0:28:99	0/1:32:98	0/0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1/1:15:45	0/1:14:52	1/1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0/1:45:99	1/0:42:99	0/1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1/1:10:35	0/0:8:30	1/1:9:32
EOF
fi

if [ ! -f data/concordance_multi_allelic.vcf ]; then
  cat > data/concordance_multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G,T	100	PASS	DP=55	GT:DP:GQ	0/1:30:99	0/1:25:95	0/1:28:98
chr1	20000	rs456	C	G,T,A	80	PASS	DP=60	GT:DP:GQ	0/2:28:99	0/2:32:98	0/2:27:97
chr1	30000	rs789	G	A,C	50	LowQual	DP=29	GT:DP:GQ	1/2:15:45	2/1:14:52	1/2:16:48
chr2	15000	rs101	T	C,A	90	PASS	DP=87	GT:DP:GQ	0/1:45:99	0/2:42:99	0/1:44:99
chr2	25000	rs102	G	T,C	30	LowQual	DP=18	GT:DP:GQ	1/1:10:35	1/1:8:30	1/1:9:32
EOF
fi

if [ ! -f data/concordance_missing_data.vcf ]; then
  cat > data/concordance_missing_data.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0/1:30:99	./.:25:95	0/1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	./.:28:99	0/0:32:98	0/0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1/1:15:45	1/1:14:52	./.:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	./.:45:99	./.:42:99	./.:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1/1:10:35	1/1:8:30	1/1:9:32
EOF
fi

if [ ! -f data/concordance_phased.vcf ]; then
  cat > data/concordance_phased.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1|1:15:45	1|1:14:52	1|1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0|1:45:99	1|0:42:99	0|1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1|1:10:35	1|1:8:30	1|1:9:32
EOF
fi

if [ ! -f data/concordance_malformed.vcf ]; then
  cat > data/concordance_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT	0/1	0/1	0/1
not_a_chromosome	20000	rs456	C	T	80	PASS	DP=60	GT	0/0	0/0	0/0
chr1	not_a_position	rs789	G	A	50	LowQual	DP=29	GT	1/1	1/1	1/1
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT	1/X	1/1	1/1
chr1	40000	rs789	G	A	50	LowQual	DP=29	GT	9/9	1/1	1/1
malformed_line_with_few_columns
EOF
fi

# Create large test data with 500 variants
if [ ! -f data/concordance_large.vcf ]; then
  cat > data/concordance_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
EOF

  # Add 500 variant lines
  for i in $(seq 1 500); do
    case $((i % 3)) in
      0) # All concordant
        echo "chr1	$((10000 + $i * 100))	rs$i	A	G	100	PASS	DP=55	GT:DP	0/1:30	0/1:25	0/1:28" >> data/concordance_large.vcf
        ;;
      1) # Discordant
        echo "chr1	$((10000 + $i * 100))	rs$i	C	T	80	PASS	DP=60	GT:DP	0/0:28	0/1:32	0/0:27" >> data/concordance_large.vcf
        ;;
      2) # Missing data
        echo "chr1	$((10000 + $i * 100))	rs$i	G	A	50	LowQual	DP=29	GT:DP	1/1:15	./.:14	1/1:16" >> data/concordance_large.vcf
        ;;
    esac
  done
fi

# Expected outputs
if [ ! -f expected/concordance_all_match.tsv ]; then
  cat > expected/concordance_all_match.tsv << EOF
CHROM	POS	ID	REF	ALT	Num_Samples	Unique_Normalized_Genotypes	Concordance_Status
chr1	10000	rs123	A	G	3	1	CONCORDANT
chr1	20000	rs456	C	T	3	1	CONCORDANT
chr1	30000	rs789	G	A	3	1	CONCORDANT
chr2	15000	rs101	T	C	3	1	CONCORDANT
chr2	25000	rs102	G	T	3	1	CONCORDANT
EOF
fi

if [ ! -f expected/concordance_some_mismatch.tsv ]; then
  cat > expected/concordance_some_mismatch.tsv << EOF
CHROM	POS	ID	REF	ALT	Num_Samples	Unique_Normalized_Genotypes	Concordance_Status
chr1	10000	rs123	A	G	3	1	CONCORDANT
chr1	20000	rs456	C	T	3	2	DISCORDANT
chr1	30000	rs789	G	A	3	2	DISCORDANT
chr2	15000	rs101	T	C	3	1	CONCORDANT
chr2	25000	rs102	G	T	3	2	DISCORDANT
EOF
fi

if [ ! -f expected/concordance_multi_allelic.tsv ]; then
  cat > expected/concordance_multi_allelic.tsv << EOF
CHROM	POS	ID	REF	ALT	Num_Samples	Unique_Normalized_Genotypes	Concordance_Status
chr1	10000	rs123	A	G,T	3	1	CONCORDANT
chr1	20000	rs456	C	G,T,A	3	1	CONCORDANT
chr1	30000	rs789	G	A,C	3	1	CONCORDANT
chr2	15000	rs101	T	C,A	3	2	DISCORDANT
chr2	25000	rs102	G	T,C	3	1	CONCORDANT
EOF
fi

if [ ! -f expected/concordance_missing_data.tsv ]; then
  cat > expected/concordance_missing_data.tsv << EOF
CHROM	POS	ID	REF	ALT	Num_Samples	Unique_Normalized_Genotypes	Concordance_Status
chr1	10000	rs123	A	G	2	1	CONCORDANT
chr1	20000	rs456	C	T	2	1	CONCORDANT
chr1	30000	rs789	G	A	2	1	CONCORDANT
chr2	15000	rs101	T	C	0	0	NO_GENOTYPES
chr2	25000	rs102	G	T	3	1	CONCORDANT
EOF
fi

if [ ! -f expected/concordance_phased.tsv ]; then
  cat > expected/concordance_phased.tsv << EOF
CHROM	POS	ID	REF	ALT	Num_Samples	Unique_Normalized_Genotypes	Concordance_Status
chr1	10000	rs123	A	G	3	1	CONCORDANT
chr1	20000	rs456	C	T	3	1	CONCORDANT
chr1	30000	rs789	G	A	3	1	CONCORDANT
chr2	15000	rs101	T	C	3	1	CONCORDANT
chr2	25000	rs102	G	T	3	1	CONCORDANT
EOF
fi

# Test 1: All samples matching
run_test 1 "All samples matching" \
  "$VCFX_EXECUTABLE < data/concordance_all_match.vcf" \
  "expected/concordance_all_match.tsv" \
  "out/concordance_all_match.tsv"

# Test 2: Some samples mismatching
run_test 2 "Some samples mismatching" \
  "$VCFX_EXECUTABLE < data/concordance_some_mismatch.vcf" \
  "expected/concordance_some_mismatch.tsv" \
  "out/concordance_some_mismatch.tsv"

# Test 3: Multi-allelic variants
run_test 3 "Multi-allelic variants" \
  "$VCFX_EXECUTABLE < data/concordance_multi_allelic.vcf" \
  "expected/concordance_multi_allelic.tsv" \
  "out/concordance_multi_allelic.tsv"

# Test 4: Missing data (./:)
run_test 4 "Missing data" \
  "$VCFX_EXECUTABLE < data/concordance_missing_data.vcf" \
  "expected/concordance_missing_data.tsv" \
  "out/concordance_missing_data.tsv"

# Test 5: Phased genotypes (|)
run_test 5 "Phased genotypes" \
  "$VCFX_EXECUTABLE < data/concordance_phased.vcf" \
  "expected/concordance_phased.tsv" \
  "out/concordance_phased.tsv"

# Test 6: Malformed input (should handle gracefully)
echo "Test 6: Malformed input (should handle gracefully)"
$VCFX_EXECUTABLE < data/concordance_malformed.vcf > out/concordance_malformed.tsv
# Just check that we got some output with expected headers
if ! grep -q "CHROM" out/concordance_malformed.tsv || ! grep -q "Concordance_Status" out/concordance_malformed.tsv; then
  echo "  Test 6 failed. Output doesn't contain expected headers."
  exit 1
fi
echo "  Test 6 passed."

# Test 7: Performance with large file
echo "Test 7: Performance with large file"
time $VCFX_EXECUTABLE < data/concordance_large.vcf > out/concordance_large.tsv
# Just check that output is non-empty and contains expected format
if ! grep -q "CHROM" out/concordance_large.tsv || ! grep -q "Concordance_Status" out/concordance_large.tsv; then
  echo "  Test 7 failed. Output doesn't contain expected headers."
  exit 1
fi

# Count the concordant/discordant/missing variants
total_variants=$(grep -v "CHROM" out/concordance_large.tsv | wc -l)
concordant_count=$(grep -c "CONCORDANT" out/concordance_large.tsv || true)
discordant_count=$(grep -c "DISCORDANT" out/concordance_large.tsv || true)
missing_count=$(grep -c "NO_GENOTYPES" out/concordance_large.tsv || true)

if [ "$total_variants" -lt 490 ]; then
  echo "  Test 7 failed. Expected at least 490 variants, but found $total_variants."
  exit 1
fi

echo "  Test 7 passed. Found $concordant_count concordant, $discordant_count discordant, and $missing_count with no genotypes."

# Test 8: Sample subset option
run_test 8 "Sample subset" \
  "$VCFX_EXECUTABLE --samples SAMPLE1,SAMPLE3 < data/concordance_some_mismatch.vcf" \
  "expected/concordance_subset.tsv" \
  "out/concordance_subset.tsv"

# Test 8: Help display
echo "Test 9: Verifying help display"
$VCFX_EXECUTABLE --help > out/concordance_help.txt
if ! grep -q "VCFX_cross_sample_concordance" out/concordance_help.txt; then
  echo "  Test 9 failed. Help message not displayed correctly."
  cat out/concordance_help.txt
  exit 1
fi
echo "  Test 9 passed."

echo "All VCFX_cross_sample_concordance tests passed!"
