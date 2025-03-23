#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_header_parser ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_header_parser/VCFX_header_parser"

# Check if executable exists
if [ ! -f "$VCFX_EXECUTABLE" ]; then
  echo "Error: $VCFX_EXECUTABLE not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Create test data if it doesn't exist
if [ ! -f data/header_parser_input.vcf ]; then
  cat > data/header_parser_input.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	80	PASS	AF=0.10	GT:DP	0/0:28	0/1:32
EOF
fi

# Create expected output
if [ ! -f expected/header_parser_output.txt ]; then
  cat > expected/header_parser_output.txt << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
EOF
fi

# Run the test
echo "Testing basic header extraction..."
$VCFX_EXECUTABLE < data/header_parser_input.vcf > out/header_parser_output.txt

# Compare results
diff -u expected/header_parser_output.txt out/header_parser_output.txt
echo "VCFX_header_parser test passed!" 