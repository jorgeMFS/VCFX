#!/usr/bin/env bash
# Remove set -e to allow continuing after errors

# Track failures
failures=0

echo "=== Testing VCFX_multiallelic_splitter ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_multiallelic_splitter/VCFX_multiallelic_splitter"

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

# Create test data: Basic multi-allelic site
if [ ! -f data/multiallelic_basic.vcf ]; then
  cat > data/multiallelic_basic.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=BaseQRankSum,Number=A,Type=Float,Description="Z-score from Wilcoxon rank sum test of alt vs. ref base qualities">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G,T	100	PASS	DP=50;AF=0.2,0.1;BaseQRankSum=1.5,2.0	GT:AD:DP:GQ:PL	0/1:10,5,0:15:20:20,0,100,40,100,100	0/2:12,0,4:16:18:18,50,100,0,100,20
EOF
fi

# Create expected output for basic test
if [ ! -f expected/multiallelic_basic_split.vcf ]; then
  cat > expected/multiallelic_basic_split.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=BaseQRankSum,Number=A,Type=Float,Description="Z-score from Wilcoxon rank sum test of alt vs. ref base qualities">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=50;AF=0.2;BaseQRankSum=1.5	GT:AD:DP:GQ:PL	0/1:10,5:15:20:20,0,100	.:.:.:.:.
chr1	10000	rs123	A	T	100	PASS	DP=50;AF=0.1;BaseQRankSum=2.0	GT:AD:DP:GQ:PL	.:.:.:.:.	0/1:12,4:16:18:18,0,100
EOF
fi

# Create test data: Complex multi-allelic sites with different INFO/FORMAT fields
if [ ! -f data/multiallelic_complex.vcf ]; then
  cat > data/multiallelic_complex.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation of allele frequency">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G,T,C	100	PASS	DP=80;AF=0.2,0.1,0.05;AC=10,5,3;MLEAF=0.18,0.11,0.06	GT:AD:DP:GQ:PL	0/1:30,10,0,0:40:30:30,0,100,40,100,100,50,120,130,150	0/2:35,0,8,0:43:25:25,50,100,0,100,20,60,110,115,140	1/3:25,12,0,7:44:35:35,0,120,50,140,160,0,130,150,0
chr2	20000	rs456	CT	C,CTT,CTTT	90	PASS	DP=70;AF=0.15,0.08,0.04;AC=8,4,2;MLEAF=0.13,0.09,0.03	GT:AD:DP:GQ:PL	0/1:25,8,0,0:33:22:22,0,90,35,90,90,40,110,115,130	0/3:30,0,0,5:35:18:18,40,90,45,95,110,0,80,90,0	2/3:22,0,6,4:32:28:28,35,85,0,80,105,0,90,0,0
chr3	30000	sv1	G	<DEL>,<INS>	80	PASS	DP=60;AF=0.12,0.06;AC=6,3;MLEAF=0.11,0.07;END=40000;SVTYPE=DEL	GT:AD:DP:GQ:PL	0/1:20,6,0:26:15:15,0,80,25,80,90	0/2:18,0,4:22:12:12,30,70,0,70,10	1/2:15,8,3:26:20:20,0,60,0,65,0
EOF
fi

# Create expected output for complex test
if [ ! -f expected/multiallelic_complex_split.vcf ]; then
  cat > expected/multiallelic_complex_split.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation of allele frequency">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=80;AF=0.2;AC=10;MLEAF=0.18	GT:AD:DP:GQ:PL	0/1:30,10:40:30:30,0,100	.:.:.:.:.	0/1:25,12:44:35:35,0,120
chr1	10000	rs123	A	T	100	PASS	DP=80;AF=0.1;AC=5;MLEAF=0.11	GT:AD:DP:GQ:PL	.:.:.:.:.	0/1:35,8:43:25:25,0,100	.:.:.:.:.
chr1	10000	rs123	A	C	100	PASS	DP=80;AF=0.05;AC=3;MLEAF=0.06	GT:AD:DP:GQ:PL	.:.:.:.:.	.:.:.:.:.	0/1:25,7:44:35:35,0,130
chr2	20000	rs456	CT	C	90	PASS	DP=70;AF=0.15;AC=8;MLEAF=0.13	GT:AD:DP:GQ:PL	0/1:25,8:33:22:22,0,90	.:.:.:.:.	.:.:.:.:.
chr2	20000	rs456	CT	CTT	90	PASS	DP=70;AF=0.08;AC=4;MLEAF=0.09	GT:AD:DP:GQ:PL	.:.:.:.:.	.:.:.:.:.	0/1:22,6:32:28:28,0,80
chr2	20000	rs456	CT	CTTT	90	PASS	DP=70;AF=0.04;AC=2;MLEAF=0.03	GT:AD:DP:GQ:PL	.:.:.:.:.	0/1:30,5:35:18:18,0,80	0/1:22,4:32:28:28,0,90
chr3	30000	sv1	G	<DEL>	80	PASS	DP=60;AF=0.12;AC=6;MLEAF=0.11;END=40000;SVTYPE=DEL	GT:AD:DP:GQ:PL	0/1:20,6:26:15:15,0,80	.:.:.:.:.	0/1:15,8:26:20:20,0,60
chr3	30000	sv1	G	<INS>	80	PASS	DP=60;AF=0.06;AC=3;MLEAF=0.07;END=40000;SVTYPE=DEL	GT:AD:DP:GQ:PL	.:.:.:.:.	0/1:18,4:22:12:12,0,70	0/1:15,3:26:20:20,0,65
EOF
fi

# Create test data: Missing values and edge cases
if [ ! -f data/multiallelic_edge.vcf ]; then
  cat > data/multiallelic_edge.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G,T	100	PASS	DP=50;AF=0.2,0.1	GT:AD:DP	./.:.:.	0/2:10,0,5:15
chr1	20000	rs456	C	G,T	80	PASS	DP=40;AF=.,0.1	GT:AD:DP	0/1:15,8,.:23	0/2:12,.,6:18
chr1	30000	.	ACTG	A,ACTGG	70	PASS	DP=35	GT:AD:DP	0/1:10,5,.:15	0/2:8,.,4:12
chr1	40000	rs789	A	G,*	60	PASS	DP=30	GT:AD:DP	0/1:10,6,.:16	0/2:9,.,3:12
chr1	50000	multi	G	<NON_REF>,T	50	PASS	DP=25	GT:AD:DP	0/2:8,.,5:13	0/1:7,4,.:11
EOF
fi

# Create expected output for edge cases test
if [ ! -f expected/multiallelic_edge_split.vcf ]; then
  cat > expected/multiallelic_edge_split.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=50;AF=0.2	GT:AD:DP	.:.:.	.:.:.
chr1	10000	rs123	A	T	100	PASS	DP=50;AF=0.1	GT:AD:DP	.:.:.	0/1:10,5:15
chr1	20000	rs456	C	G	80	PASS	DP=40;AF=.	GT:AD:DP	0/1:15,8:23	.:.:.
chr1	20000	rs456	C	T	80	PASS	DP=40;AF=0.1	GT:AD:DP	.:.:.	0/1:12,6:18
chr1	30000	.	ACTG	A	70	PASS	DP=35	GT:AD:DP	0/1:10,5:15	.:.:.
chr1	30000	.	ACTG	ACTGG	70	PASS	DP=35	GT:AD:DP	.:.:.	0/1:8,4:12
chr1	40000	rs789	A	G	60	PASS	DP=30	GT:AD:DP	0/1:10,6:16	.:.:.
chr1	40000	rs789	A	*	60	PASS	DP=30	GT:AD:DP	.:.:.	0/1:9,3:12
chr1	50000	multi	G	<NON_REF>	50	PASS	DP=25	GT:AD:DP	.:.:.	0/1:7,4:11
chr1	50000	multi	G	T	50	PASS	DP=25	GT:AD:DP	0/1:8,5:13	.:.:.
EOF
fi

# Create test data: Phased genotypes
if [ ! -f data/multiallelic_phased.vcf ]; then
  cat > data/multiallelic_phased.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G,T	100	PASS	DP=50;AF=0.2,0.1	GT:AD:DP	0|1:10,5,0:15	0|2:12,0,4:16
chr1	20000	rs456	C	G,T	80	PASS	DP=40;AF=0.15,0.08	GT:AD:DP	1|0:15,8,0:23	0|2:12,0,5:17
EOF
fi

# Create expected output for phased test
if [ ! -f expected/multiallelic_phased_split.vcf ]; then
  cat > expected/multiallelic_phased_split.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	DP=50;AF=0.2	GT:AD:DP	0/1:10,5:15	.:.:.
chr1	10000	rs123	A	T	100	PASS	DP=50;AF=0.1	GT:AD:DP	.:.:.	0/1:12,4:16
chr1	20000	rs456	C	G	80	PASS	DP=40;AF=0.15	GT:AD:DP	0/1:15,8:23	.:.:.
chr1	20000	rs456	C	T	80	PASS	DP=40;AF=0.08	GT:AD:DP	.:.:.	0/1:12,5:17
EOF
fi

# Run the tests
# Test 1: Basic multi-allelic site with simple genotypes
run_test 1 "Basic multi-allelic site splitting" \
  "$VCFX_EXECUTABLE < data/multiallelic_basic.vcf" \
  "expected/multiallelic_basic_split.vcf" \
  "out/multiallelic_basic_split.vcf"

# Test 2: Complex multi-allelic sites with multiple ALTs and samples
run_test 2 "Complex multi-allelic sites with multiple ALTs" \
  "$VCFX_EXECUTABLE < data/multiallelic_complex.vcf" \
  "expected/multiallelic_complex_split.vcf" \
  "out/multiallelic_complex_split.vcf"

# Test 3: Missing values and edge cases
run_test 3 "Edge cases with missing values" \
  "$VCFX_EXECUTABLE < data/multiallelic_edge.vcf" \
  "expected/multiallelic_edge_split.vcf" \
  "out/multiallelic_edge_split.vcf"

# Test 4: Phased genotypes
run_test 4 "Phased genotypes" \
  "$VCFX_EXECUTABLE < data/multiallelic_phased.vcf" \
  "expected/multiallelic_phased_split.vcf" \
  "out/multiallelic_phased_split.vcf"

# Test 5: Help message
echo "Test 5: Help message"
$VCFX_EXECUTABLE --help > out/multiallelic_help.txt
if ! grep -q "VCFX_multiallelic_splitter" out/multiallelic_help.txt; then
  echo "  Test 5 failed: Help message doesn't contain expected content"
  ((failures++))
else
  echo "  Test 5 passed."
fi

# Report results
if [ $failures -eq 0 ]; then
  echo "All VCFX_multiallelic_splitter tests passed!"
  exit 0
else
  echo "$failures VCFX_multiallelic_splitter tests failed."
  exit 1
fi 