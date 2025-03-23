#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_field_extractor ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_field_extractor/VCFX_field_extractor"

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
if [ ! -f data/field_extractor_input.vcf ]; then
  cat > data/field_extractor_input.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25;TYPE=SNP	GT:DP:GQ	0/1:30:99	0/0:25:95
chr1	20000	rs456	C	T	80	PASS	DP=60;AF=0.10;TYPE=SNP	GT:DP:GQ	0/0:28:99	0/1:32:98
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05;TYPE=SNP	GT:DP:GQ	0/1:15:45	0/0:14:52
chr1	40000	.	T	C	90	PASS	DP=87;AF=0.35;TYPE=SNP	GT:DP:GQ	0/1:45:99	1/1:42:99
chr1	50000	rs102	G	T	30	LowQual	DP=18;AF=0.02;TYPE=SNP	GT:DP:GQ	0/0:10:35	0/1:8:30
EOF
fi

# Create malformed VCF file
if [ ! -f data/field_extractor_malformed.vcf ]; then
  cat > data/field_extractor_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP:GQ	0/1:30	0/0:25
chr1	invalid_pos	rs456	C	T	80	PASS	DP=60;AF=INVALID	GT:DP:GQ	0/0:28:99	0/1:32:98
malformed_line_with_fewer_columns	20000	rs789
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05	INVALID_FORMAT	0/1:15:45	0/0:14:52
EOF
fi

# Create a large VCF file for performance testing
if [ ! -f data/field_extractor_large.vcf ]; then
  cat > data/field_extractor_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF

  # Add 500 variant lines
  for i in $(seq 1 500); do
    echo "chr1	$((10000 + $i * 100))	rs$i	A	G	100	PASS	DP=$((50 + $i));AF=$((i % 10 + 1))0.01;TYPE=SNP	GT:DP:GQ	0/1:$((30 + $i % 20)):99	0/0:$((25 + $i % 15)):95" >> data/field_extractor_large.vcf
  done
fi

# Create a file with headers only (no variants)
if [ ! -f data/field_extractor_headers_only.vcf ]; then
  cat > data/field_extractor_headers_only.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF
fi

# Test 1: Basic field extraction
if [ ! -f expected/field_basic.tsv ]; then
  cat > expected/field_basic.tsv << EOF
CHROM	POS	ID	REF	ALT
chr1	10000	rs123	A	G
chr1	20000	rs456	C	T
chr1	30000	rs789	G	A
chr1	40000	.	T	C
chr1	50000	rs102	G	T
EOF
fi

run_test 1 "Basic field extraction (standard fields)" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,ID,REF,ALT\" < data/field_extractor_input.vcf" \
  "expected/field_basic.tsv" \
  "out/field_basic.tsv"

# Test 2: Extracting INFO fields
if [ ! -f expected/field_info.tsv ]; then
  cat > expected/field_info.tsv << EOF
CHROM	POS	DP	AF	TYPE
chr1	10000	55	0.25	SNP
chr1	20000	60	0.10	SNP
chr1	30000	29	0.05	SNP
chr1	40000	87	0.35	SNP
chr1	50000	18	0.02	SNP
EOF
fi

run_test 2 "INFO field extraction" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,DP,AF,TYPE\" < data/field_extractor_input.vcf" \
  "expected/field_info.tsv" \
  "out/field_info.tsv"

# Test 3: Extracting sample fields
if [ ! -f expected/field_samples.tsv ]; then
  cat > expected/field_samples.tsv << EOF
CHROM	POS	SAMPLE1:GT	SAMPLE1:DP	SAMPLE2:GT	SAMPLE2:DP
chr1	10000	0/1	30	0/0	25
chr1	20000	0/0	28	0/1	32
chr1	30000	0/1	15	0/0	14
chr1	40000	0/1	45	1/1	42
chr1	50000	0/0	10	0/1	8
EOF
fi

run_test 3 "Sample field extraction" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,SAMPLE1:GT,SAMPLE1:DP,SAMPLE2:GT,SAMPLE2:DP\" < data/field_extractor_input.vcf" \
  "expected/field_samples.tsv" \
  "out/field_samples.tsv"

# Test 4: Mixed fields extraction
if [ ! -f expected/field_mixed.tsv ]; then
  cat > expected/field_mixed.tsv << EOF
CHROM	POS	ID	DP	AF	SAMPLE1:GT	SAMPLE2:GT
chr1	10000	rs123	55	0.25	0/1	0/0
chr1	20000	rs456	60	0.10	0/0	0/1
chr1	30000	rs789	29	0.05	0/1	0/0
chr1	40000	.	87	0.35	0/1	1/1
chr1	50000	rs102	18	0.02	0/0	0/1
EOF
fi

run_test 4 "Mixed fields extraction (standard, INFO, and sample fields)" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,ID,DP,AF,SAMPLE1:GT,SAMPLE2:GT\" < data/field_extractor_input.vcf" \
  "expected/field_mixed.tsv" \
  "out/field_mixed.tsv"

# Test 5: Sample indexing with S1, S2 format
if [ ! -f expected/field_sample_index.tsv ]; then
  cat > expected/field_sample_index.tsv << EOF
CHROM	POS	S1:GT	S1:DP	S2:GT	S2:DP
chr1	10000	0/1	30	0/0	25
chr1	20000	0/0	28	0/1	32
chr1	30000	0/1	15	0/0	14
chr1	40000	0/1	45	1/1	42
chr1	50000	0/0	10	0/1	8
EOF
fi

run_test 5 "Sample fields using S1, S2 notation" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,S1:GT,S1:DP,S2:GT,S2:DP\" < data/field_extractor_input.vcf" \
  "expected/field_sample_index.tsv" \
  "out/field_sample_index.tsv"

# Test 6: Non-existent fields
if [ ! -f expected/field_nonexistent.tsv ]; then
  cat > expected/field_nonexistent.tsv << EOF
CHROM	POS	NONEXISTENT_FIELD	NONEXISTENT_INFO	NONEXISTENT_SAMPLE:GT
chr1	10000	.	.	.
chr1	20000	.	.	.
chr1	30000	.	.	.
chr1	40000	.	.	.
chr1	50000	.	.	.
EOF
fi

run_test 6 "Non-existent fields (should return '.' placeholders)" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,NONEXISTENT_FIELD,NONEXISTENT_INFO,NONEXISTENT_SAMPLE:GT\" < data/field_extractor_input.vcf" \
  "expected/field_nonexistent.tsv" \
  "out/field_nonexistent.tsv"

# Test 7: Malformed VCF file
if [ ! -f expected/field_malformed.tsv ]; then
  cat > expected/field_malformed.tsv << EOF
CHROM	POS	ID	INFO	SAMPLE1:GT
chr1	10000	rs123	DP=55;AF=0.25	0/1
chr1	invalid_pos	rs456	DP=60;AF=INVALID	0/0
malformed_line_with_fewer_columns	20000	rs789	.	.
chr1	30000	rs789	DP=29;AF=0.05	.
EOF
fi

run_test 7 "Malformed VCF file (should handle gracefully)" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,ID,INFO,SAMPLE1:GT\" < data/field_extractor_malformed.vcf" \
  "expected/field_malformed.tsv" \
  "out/field_malformed.tsv"

# Test 8: Empty/headers-only VCF file
if [ ! -f expected/field_headers_only.tsv ]; then
  cat > expected/field_headers_only.tsv << EOF
CHROM	POS	ID	REF	ALT
EOF
fi

run_test 8 "Headers-only VCF file (should produce header row only)" \
  "$VCFX_EXECUTABLE --fields \"CHROM,POS,ID,REF,ALT\" < data/field_extractor_headers_only.vcf" \
  "expected/field_headers_only.tsv" \
  "out/field_headers_only.tsv"

# Test 9: Performance with large VCF file
if [ ! -f expected/field_large.tsv ]; then
  # Create only the first few lines of expected output to verify
  cat > expected/field_large_sample.tsv << EOF
CHROM	POS	ID	DP	AF	SAMPLE1:GT	SAMPLE2:GT
chr1	10100	rs1	51	20.01	0/1	0/0
chr1	10200	rs2	52	30.01	0/1	0/0
chr1	10300	rs3	53	40.01	0/1	0/0
EOF
fi

echo "Test 9: Performance with large VCF file"
time $VCFX_EXECUTABLE --fields "CHROM,POS,ID,DP,AF,SAMPLE1:GT,SAMPLE2:GT" < data/field_extractor_large.vcf > out/field_large.tsv
# Verify first few lines only for the large file
head -n 4 out/field_large.tsv > out/field_large_sample.tsv
diff -u expected/field_large_sample.tsv out/field_large_sample.tsv || {
  echo "  Test 9 failed. Actual output (first few lines):"
  cat out/field_large_sample.tsv
  exit 1
}
echo "  Test 9 passed."

# Test 10: Help display
echo "Test 10: Verifying help display"
$VCFX_EXECUTABLE --help > out/field_extractor_help.txt
if ! grep -q "VCFX_field_extractor" out/field_extractor_help.txt; then
  echo "  Test 10 failed. Help message not displayed correctly."
  cat out/field_extractor_help.txt
  exit 1
fi
echo "  Test 10 passed."

echo "All VCFX_field_extractor tests passed!" 