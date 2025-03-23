#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_sample_extractor ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_sample_extractor/VCFX_sample_extractor"

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

# Create test data if it doesn't exist
if [ ! -f data/sample_extractor_input.vcf ]; then
  cat > data/sample_extractor_input.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30	0/0:25	1/1:20	0/1:22
chr1	20000	rs456	C	T	80	PASS	DP=60;AF=0.10	GT:DP	0/0:28	0/1:32	0/0:27	1/1:30
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05	GT:DP	0/1:15	0/0:14	0/1:12	0/0:13
EOF
fi

if [ ! -f expected/sample_extractor_single.vcf ]; then
  cat > expected/sample_extractor_single.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/1:30
chr1	20000	rs456	C	T	80	PASS	DP=60;AF=0.10	GT:DP	0/0:28
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05	GT:DP	0/1:15
EOF
fi

if [ ! -f expected/sample_extractor_multiple.vcf ]; then
  cat > expected/sample_extractor_multiple.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE2	SAMPLE4
chr1	10000	rs123	A	G	100	PASS	DP=55;AF=0.25	GT:DP	0/0:25	0/1:22
chr1	20000	rs456	C	T	80	PASS	DP=60;AF=0.10	GT:DP	0/1:32	1/1:30
chr1	30000	rs789	G	A	50	LowQual	DP=29;AF=0.05	GT:DP	0/0:14	0/0:13
EOF
fi

# Test 1: Extract a single sample
echo "Test 1: Extracting a single sample"
$VCFX_EXECUTABLE --samples "SAMPLE1" < data/sample_extractor_input.vcf > out/sample_extractor_single.vcf
diff -u expected/sample_extractor_single.vcf out/sample_extractor_single.vcf
echo "  Test 1 passed."

# Test 2: Extract multiple samples
echo "Test 2: Extracting multiple samples"
$VCFX_EXECUTABLE --samples "SAMPLE2,SAMPLE4" < data/sample_extractor_input.vcf > out/sample_extractor_multiple.vcf
diff -u expected/sample_extractor_multiple.vcf out/sample_extractor_multiple.vcf
echo "  Test 2 passed."

# Test 3: Extract with space-delimited sample names
echo "Test 3: Extracting with space-delimited sample names"
$VCFX_EXECUTABLE --samples "SAMPLE2 SAMPLE4" < data/sample_extractor_input.vcf > out/sample_extractor_multiple_space.vcf
diff -u expected/sample_extractor_multiple.vcf out/sample_extractor_multiple_space.vcf
echo "  Test 3 passed."

echo "All VCFX_sample_extractor tests passed!" 