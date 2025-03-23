#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_variant_classifier ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_variant_classifier/VCFX_variant_classifier"

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
if [ ! -f data/classifier_mixed.vcf ]; then
  cat > data/classifier_mixed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	snp1	A	G	100	PASS	DP=55;AF=0.25	GT	0/1
chr1	20000	indel1	CA	T	80	PASS	DP=60;AF=0.10	GT	0/1
chr1	30000	indel2	G	GACGT	50	LowQual	DP=29;AF=0.05	GT	0/1
chr1	40000	mnv1	ATG	GTC	90	PASS	DP=87;AF=0.35	GT	0/1
chr1	50000	sv1	G	<DEL>	30	LowQual	SVTYPE=DEL;END=51000	GT	0/1
chr2	15000	sv2	T	T[chr3:12345[	70	PASS	DP=40;AF=0.15	GT	0/1
chr2	25000	sv3	AACTGAATGACTGACTGACTGACTGACTAGCTAGCTAGCTAGCTAGCTA	A	60	PASS	DP=30;AF=0.20	GT	0/1
chr2	35000	sv4	A	AACTGAATGACTGACTGACTGACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA	50	PASS	DP=25;AF=0.18	GT	0/1
EOF
fi

# Create multi-allelic test data
if [ ! -f data/classifier_multi_allelic.vcf ]; then
  cat > data/classifier_multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	multi1	A	G,T	100	PASS	DP=55;AF=0.25	GT	0/1
chr1	20000	multi2	CAT	C,CATGGT	80	PASS	DP=60;AF=0.10	GT	1/2
chr1	30000	multi3	G	GACGT,<DEL>,GG	50	LowQual	DP=29;AF=0.05	GT	0/2
chr1	40000	multi4	ATG	GTC,A,ATGC	90	PASS	DP=87;AF=0.35	GT	0/3
EOF
fi

# Create edge case test data
if [ ! -f data/classifier_edge_cases.vcf ]; then
  cat > data/classifier_edge_cases.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	edge1	A	.	100	PASS	DP=55	GT	0/0
chr1	20000	edge2	G	G	80	PASS	DP=60	GT	0/0
chr1	30000	edge3	T	<NON_REF>	50	LowQual	DP=29	GT	0/1
chr1	40000	edge4	C	*	90	PASS	DP=87	GT	0/1
chr1	50000	edge5	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	A	60	PASS	DP=30	GT	0/1
chr1	60000	edge6	A	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	70	PASS	DP=40	GT	0/1
EOF
fi

# Create malformed test data
if [ ! -f data/classifier_malformed.vcf ]; then
  cat > data/classifier_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	mal1	A	G	100	PASS	.
malformed_line_with_few_columns
chr1	not_a_position	mal2	C	T	80	PASS	.
chr1	30000	mal3	G		50	LowQual	.
chr1	40000	mal4		A	90	PASS	.
chr1	50000	mal5	G	T,	30	LowQual	.
EOF
fi

# Create large test data for performance testing
if [ ! -f data/classifier_large.vcf ]; then
  cat > data/classifier_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
EOF
  
  # Add 500 variant lines with different types
  for i in $(seq 1 500); do
    case $((i % 4)) in
      0) # SNP
        echo "chr1	$((10000 + $i * 100))	var$i	A	G	100	PASS	DP=$((50 + $i))	GT	0/1" >> data/classifier_large.vcf
        ;;
      1) # INDEL
        echo "chr1	$((10000 + $i * 100))	var$i	A	AT	100	PASS	DP=$((50 + $i))	GT	0/1" >> data/classifier_large.vcf
        ;;
      2) # MNV
        echo "chr1	$((10000 + $i * 100))	var$i	AT	GC	100	PASS	DP=$((50 + $i))	GT	0/1" >> data/classifier_large.vcf
        ;;
      3) # STRUCTURAL
        echo "chr1	$((10000 + $i * 100))	var$i	A	<DEL>	100	PASS	DP=$((50 + $i))	GT	0/1" >> data/classifier_large.vcf
        ;;
    esac
  done
fi

# Expected outputs
if [ ! -f expected/classifier_mixed.tsv ]; then
  cat > expected/classifier_mixed.tsv << EOF
CHROM	POS	ID	REF	ALT	Classification
chr1	10000	snp1	A	G	SNP
chr1	20000	indel1	CA	T	INDEL
chr1	30000	indel2	G	GACGT	INDEL
chr1	40000	mnv1	ATG	GTC	MNV
chr1	50000	sv1	G	<DEL>	STRUCTURAL
chr2	15000	sv2	T	T[chr3:12345[	STRUCTURAL
chr2	25000	sv3	AACTGAATGACTGACTGACTGACTAGCTAGCTAGCTAGCTAGCTA	A	STRUCTURAL
chr2	35000	sv4	A	AACTGAATGACTGACTGACTGACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA	STRUCTURAL
EOF
fi

if [ ! -f expected/classifier_multi_allelic.tsv ]; then
  cat > expected/classifier_multi_allelic.tsv << EOF
CHROM	POS	ID	REF	ALT	Classification
chr1	10000	multi1	A	G,T	SNP
chr1	20000	multi2	CAT	C,CATGGT	INDEL
chr1	30000	multi3	G	GACGT,<DEL>,GG	STRUCTURAL
chr1	40000	multi4	ATG	GTC,A,ATGC	INDEL
EOF
fi

if [ ! -f expected/classifier_edge_cases.tsv ]; then
  cat > expected/classifier_edge_cases.tsv << EOF
CHROM	POS	ID	REF	ALT	Classification
chr1	10000	edge1	A	.	UNKNOWN
chr1	20000	edge2	G	G	UNKNOWN
chr1	30000	edge3	T	<NON_REF>	STRUCTURAL
chr1	40000	edge4	C	*	UNKNOWN
chr1	50000	edge5	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	A	STRUCTURAL
chr1	60000	edge6	A	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	STRUCTURAL
EOF
fi

if [ ! -f expected/classifier_malformed.tsv ]; then
  cat > expected/classifier_malformed.tsv << EOF
CHROM	POS	ID	REF	ALT	Classification
chr1	10000	mal1	A	G	SNP
EOF
fi

if [ ! -f expected/classifier_mixed_append.vcf ]; then
  cat > expected/classifier_mixed_append.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	snp1	A	G	100	PASS	DP=55;AF=0.25;VCF_CLASS=SNP	GT	0/1
chr1	20000	indel1	CA	T	80	PASS	DP=60;AF=0.10;VCF_CLASS=INDEL	GT	0/1
chr1	30000	indel2	G	GACGT	50	LowQual	DP=29;AF=0.05;VCF_CLASS=INDEL	GT	0/1
chr1	40000	mnv1	ATG	GTC	90	PASS	DP=87;AF=0.35;VCF_CLASS=MNV	GT	0/1
chr1	50000	sv1	G	<DEL>	30	LowQual	SVTYPE=DEL;END=51000;VCF_CLASS=STRUCTURAL	GT	0/1
chr2	15000	sv2	T	T[chr3:12345[	70	PASS	DP=40;AF=0.15;VCF_CLASS=STRUCTURAL	GT	0/1
chr2	25000	sv3	AACTGAATGACTGACTGACTGACTAGCTAGCTAGCTAGCTAGCTA	A	60	PASS	DP=30;AF=0.20;VCF_CLASS=STRUCTURAL	GT	0/1
chr2	35000	sv4	A	AACTGAATGACTGACTGACTGACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA	50	PASS	DP=25;AF=0.18;VCF_CLASS=STRUCTURAL	GT	0/1
EOF
fi

# Test 1: Basic classification to TSV
run_test 1 "Basic classification to TSV" \
  "$VCFX_EXECUTABLE < data/classifier_mixed.vcf" \
  "expected/classifier_mixed.tsv" \
  "out/classifier_mixed.tsv"

# Test 2: Multi-allelic variants to TSV
run_test 2 "Multi-allelic variants to TSV" \
  "$VCFX_EXECUTABLE < data/classifier_multi_allelic.vcf" \
  "expected/classifier_multi_allelic.tsv" \
  "out/classifier_multi_allelic.tsv"

# Test 3: Edge cases to TSV
run_test 3 "Edge cases to TSV" \
  "$VCFX_EXECUTABLE < data/classifier_edge_cases.vcf" \
  "expected/classifier_edge_cases.tsv" \
  "out/classifier_edge_cases.tsv"

# Test 4: Malformed input to TSV (should handle gracefully)
run_test 4 "Malformed input to TSV" \
  "$VCFX_EXECUTABLE < data/classifier_malformed.vcf" \
  "expected/classifier_malformed.tsv" \
  "out/classifier_malformed.tsv"

# Test 5: Append classification to INFO field
run_test 5 "Append classification to INFO field" \
  "$VCFX_EXECUTABLE --append-info < data/classifier_mixed.vcf" \
  "expected/classifier_mixed_append.vcf" \
  "out/classifier_mixed_append.vcf"

# Test 6: Performance with large file (TSV mode)
echo "Test 6: Performance with large file (TSV mode)"
time $VCFX_EXECUTABLE < data/classifier_large.vcf > out/classifier_large.tsv
# Verify first few lines
head -n 4 out/classifier_large.tsv > out/classifier_large_sample.tsv
if ! grep -q "CHROM" out/classifier_large_sample.tsv || ! grep -q "Classification" out/classifier_large_sample.tsv; then
  echo "  Test 6 failed. Output doesn't contain expected headers."
  exit 1
fi

# Count the number of each variant type
total_variants=$(grep -v "CHROM" out/classifier_large.tsv | wc -l)
snp_count=$(grep -v "CHROM" out/classifier_large.tsv | grep -c "SNP" || true)
indel_count=$(grep -v "CHROM" out/classifier_large.tsv | grep -c "INDEL" || true)
mnv_count=$(grep -v "CHROM" out/classifier_large.tsv | grep -c "MNV" || true)
sv_count=$(grep -v "CHROM" out/classifier_large.tsv | grep -c "STRUCTURAL" || true)

if [ "$total_variants" -lt 490 ]; then
  echo "  Test 6 failed. Expected at least 490 classified variants, but found $total_variants."
  exit 1
fi

if [ "$snp_count" -lt 100 ] || [ "$indel_count" -lt 100 ] || [ "$mnv_count" -lt 100 ] || [ "$sv_count" -lt 100 ]; then
  echo "  Test 6 failed. Not all variant types were classified correctly."
  echo "  SNPs: $snp_count, INDELs: $indel_count, MNVs: $mnv_count, STRUCTURAL: $sv_count"
  exit 1
fi
echo "  Test 6 passed."

# Test 7: Performance with large file (append-info mode)
echo "Test 7: Performance with large file (append-info mode)"
time $VCFX_EXECUTABLE --append-info < data/classifier_large.vcf > out/classifier_large_append.vcf
# Verify first few lines
head -n 15 out/classifier_large_append.vcf > out/classifier_large_append_sample.vcf
if ! grep -q "#CHROM" out/classifier_large_append_sample.vcf; then
  echo "  Test 7 failed. Output doesn't contain VCF header."
  exit 1
fi

# Check for VCF_CLASS in the output
if ! grep -q "VCF_CLASS=SNP" out/classifier_large_append.vcf && \
   ! grep -q "VCF_CLASS=INDEL" out/classifier_large_append.vcf && \
   ! grep -q "VCF_CLASS=MNV" out/classifier_large_append.vcf && \
   ! grep -q "VCF_CLASS=STRUCTURAL" out/classifier_large_append.vcf; then
  echo "  Test 7 failed. Output doesn't contain VCF_CLASS annotations."
  exit 1
fi
echo "  Test 7 passed."

# Test 8: Help display
echo "Test 8: Verifying help display"
$VCFX_EXECUTABLE --help > out/classifier_help.txt
if ! grep -q "VCFX_variant_classifier" out/classifier_help.txt || ! grep -q "append-info" out/classifier_help.txt; then
  echo "  Test 8 failed. Help message not displayed correctly."
  cat out/classifier_help.txt
  exit 1
fi
echo "  Test 8 passed."

echo "All VCFX_variant_classifier tests passed!" 