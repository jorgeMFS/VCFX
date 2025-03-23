#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_format_converter ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_format_converter/VCFX_format_converter"

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
if [ ! -f data/format_converter_input.vcf ]; then
  cat > data/format_converter_input.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	CA	T	80	PASS	DP=60;AF=0.10	GT:DP	0/0:28	0/1:32
chr1	30000	rs789,rs999	G	A,C	50	LowQual	DP=29;AF=0.05	GT:DP	0/1:15	0/0:14
chr1	40000	.	T	C	90	PASS	DP=87;AF=0.35	GT:DP	0/1:45	1/1:42
chr1	50000	rs102	GACT	T	30	LowQual	DP=18;AF=0.02	GT:DP	0/0:10	0/1:8
EOF
fi

# Create test data for complex variants
if [ ! -f data/format_converter_complex.vcf ]; then
  cat > data/format_converter_complex.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	del1	AGCTAGCT	A	100	PASS	SVTYPE=DEL;END=10008	GT	0/1
chr2	20000	ins1	T	TAACCTTGGCC	120	PASS	SVTYPE=INS	GT	1/1
chr3	30000	dup1	G	<DUP>	80	PASS	SVTYPE=DUP;END=35000	GT	0/1
chr4	40000	inv1	T	<INV>	60	PASS	SVTYPE=INV;END=45000	GT	0/1
chrX	15000	cnv1	C	<CN0>	70	PASS	SVTYPE=CNV;END=20000	GT	1/1
EOF
fi

# Create test data with special characters for CSV testing
if [ ! -f data/format_converter_special.vcf ]; then
  cat > data/format_converter_special.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	id1	A	G	100	PASS	Note="This, has, commas"	GT	0/1
chr2	20000	id2	C	T	90	PASS	Text="Has ""quoted"" text"	GT	1/1
chr3	30000	"id,3"	G	A	80	PASS	Key="Value with ""double"" and 'single' quotes"	GT	0/1
chr4	40000	id4	T	C	70	LowQual	Descr="New
line in field"	GT	0/0
EOF
fi

# Create malformed VCF file
if [ ! -f data/format_converter_malformed.vcf ]; then
  cat > data/format_converter_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25
malformed_without_enough_columns
chr3	NOT_A_NUMBER	rs789	G	A	50	LowQual	DP=29;AF=0.05
chr4	40000	.	T	C	INVALID	PASS	DP=87;AF=0.35
EOF
fi

# Create a large VCF file for performance testing
if [ ! -f data/format_converter_large.vcf ]; then
  cat > data/format_converter_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF

  # Add 1000 variant lines
  for i in $(seq 1 1000); do
    chrom="chr$((1 + i % 22))"
    pos=$((10000 + $i * 100))
    echo "$chrom	$pos	rs$i	A	G	100	PASS	DP=$((50 + $i));AF=0.25	GT:DP	0/1:30	0/0:25" >> data/format_converter_large.vcf
  done
fi

# Create a VCF file with only headers (no variants)
if [ ! -f data/format_converter_headers_only.vcf ]; then
  cat > data/format_converter_headers_only.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF
fi

# Test 1: Convert to BED
if [ ! -f expected/format_converter_bed.bed ]; then
  cat > expected/format_converter_bed.bed << EOF
chr1	9999	10000	rs123
chr1	19999	20001	rs456
chr1	29999	30000	rs789,rs999
chr1	39999	40000	.
chr1	49999	50003	rs102
EOF
fi

run_test 1 "Converting VCF to BED format" \
  "$VCFX_EXECUTABLE --to-bed < data/format_converter_input.vcf" \
  "expected/format_converter_bed.bed" \
  "out/format_converter_bed.bed"

# Test 2: Convert to CSV
if [ ! -f expected/format_converter_csv.csv ]; then
  cat > expected/format_converter_csv.csv << EOF
CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE1,SAMPLE2
chr1,10000,rs123,A,G,100,PASS,DP=55;AF=0.25,GT:DP,0/1:30,0/0:25
chr1,20000,rs456,CA,T,80,PASS,DP=60;AF=0.10,GT:DP,0/0:28,0/1:32
chr1,30000,"rs789,rs999",G,"A,C",50,LowQual,DP=29;AF=0.05,GT:DP,0/1:15,0/0:14
chr1,40000,.,T,C,90,PASS,DP=87;AF=0.35,GT:DP,0/1:45,1/1:42
chr1,50000,rs102,GACT,T,30,LowQual,DP=18;AF=0.02,GT:DP,0/0:10,0/1:8
EOF
fi

run_test 2 "Converting VCF to CSV format" \
  "$VCFX_EXECUTABLE --to-csv < data/format_converter_input.vcf" \
  "expected/format_converter_csv.csv" \
  "out/format_converter_csv.csv"

# Test 3: Complex variants to BED
if [ ! -f expected/format_converter_complex_bed.bed ]; then
  cat > expected/format_converter_complex_bed.bed << EOF
chr1	9999	10007	del1
chr2	19999	20000	ins1
chr3	29999	30000	dup1
chr4	39999	40000	inv1
chrX	14999	15000	cnv1
EOF
fi

run_test 3 "Converting complex variants to BED format" \
  "$VCFX_EXECUTABLE --to-bed < data/format_converter_complex.vcf" \
  "expected/format_converter_complex_bed.bed" \
  "out/format_converter_complex_bed.bed"

# Test 4: Special characters to CSV
if [ ! -f expected/format_converter_special_csv.csv ]; then
  cat > expected/format_converter_special_csv.csv << EOF
CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE1
chr1,10000,id1,A,G,100,PASS,"Note=""This, has, commas""",GT,0/1
chr2,20000,id2,C,T,90,PASS,"Text=""Has """"quoted"""" text""",GT,1/1
chr3,30000,"""id,3""",G,A,80,PASS,"Key=""Value with """"double"""" and 'single' quotes""",GT,0/1
chr4,40000,id4,T,C,70,LowQual,"Descr=""New
line in field""",GT,0/0
EOF
fi

run_test 4 "Converting VCF with special characters to CSV format" \
  "$VCFX_EXECUTABLE --to-csv < data/format_converter_special.vcf" \
  "expected/format_converter_special_csv.csv" \
  "out/format_converter_special_csv.csv"

# Test 5: Malformed VCF to BED (should handle gracefully)
if [ ! -f expected/format_converter_malformed_bed.bed ]; then
  cat > expected/format_converter_malformed_bed.bed << EOF
chr1	9999	10000	rs123
chr4	39999	40000	.
EOF
fi

run_test 5 "Converting malformed VCF to BED format (should handle gracefully)" \
  "$VCFX_EXECUTABLE --to-bed < data/format_converter_malformed.vcf" \
  "expected/format_converter_malformed_bed.bed" \
  "out/format_converter_malformed_bed.bed"

# Test 6: Headers-only VCF (should produce empty output)
if [ ! -f expected/format_converter_headers_only_bed.bed ]; then
  # This should be an empty file
  touch expected/format_converter_headers_only_bed.bed
fi

run_test 6 "Converting headers-only VCF to BED format (should produce empty output)" \
  "$VCFX_EXECUTABLE --to-bed < data/format_converter_headers_only.vcf" \
  "expected/format_converter_headers_only_bed.bed" \
  "out/format_converter_headers_only_bed.bed"

# Test 7: No valid output format specified (should fail)
echo "Test 7: No valid output format specified (should fail)"
if $VCFX_EXECUTABLE < data/format_converter_input.vcf > out/format_converter_no_format.txt 2>&1; then
  echo "  Error: Expected failure when no format specified but got success"
  exit 1
else
  echo "  Test 7 passed (expected failure when no format specified)"
fi

# Test 8: Performance with large VCF file
if [ ! -f expected/format_converter_large_bed_sample.bed ]; then
  cat > expected/format_converter_large_bed_sample.bed << EOF
chr1	9999	10000	rs1
chr2	10099	10100	rs2
chr3	10199	10200	rs3
EOF
fi

echo "Test 8: Performance with large VCF file"
time $VCFX_EXECUTABLE --to-bed < data/format_converter_large.vcf > out/format_converter_large.bed
# Verify first few lines only for the large file
head -n 3 out/format_converter_large.bed > out/format_converter_large_bed_sample.bed
diff -u expected/format_converter_large_bed_sample.bed out/format_converter_large_bed_sample.bed || {
  echo "  Test 8 failed. Actual output (first few lines):"
  cat out/format_converter_large_bed_sample.bed
  exit 1
}
echo "  Test 8 passed."

# Test 9: Help display
echo "Test 9: Verifying help display"
$VCFX_EXECUTABLE --help > out/format_converter_help.txt
if ! grep -q "VCFX_format_converter" out/format_converter_help.txt; then
  echo "  Test 9 failed. Help message not displayed correctly."
  cat out/format_converter_help.txt
  exit 1
fi
echo "  Test 9 passed."

echo "All VCFX_format_converter tests passed!" 