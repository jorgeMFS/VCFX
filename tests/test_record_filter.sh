#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_record_filter ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_record_filter/VCFX_record_filter"

# Check if executable exists
if [ ! -f "$VCFX_EXECUTABLE" ]; then
  echo "Error: $VCFX_EXECUTABLE not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Create test data if it doesn't exist
if [ ! -f data/record_filter_input.vcf ]; then
  cat > data/record_filter_input.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25;DP=30	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	80	PASS	AF=0.10;DP=60	GT:DP	0/0:28	0/1:32
chr1	30000	rs789	G	A	50	LowQual	AF=0.05;DP=15	GT:DP	0/1:15	0/0:14
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
chr1	50000	rs102	G	T	30	LowQual	AF=0.02;DP=10	GT:DP	0/0:10	0/1:8
EOF
fi

# Test 1: Filter by QUAL
if [ ! -f expected/record_filter_qual.vcf ]; then
  cat > expected/record_filter_qual.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25;DP=30	GT:DP	0/1:30	0/0:25
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
EOF
fi

echo "Test 1: Filtering by QUAL>=90"
$VCFX_EXECUTABLE --filter "QUAL>=90" < data/record_filter_input.vcf > out/record_filter_qual.vcf
diff -u expected/record_filter_qual.vcf out/record_filter_qual.vcf
echo "  Test 1 passed."

# Test 2: Filter by INFO field
if [ ! -f expected/record_filter_af.vcf ]; then
  cat > expected/record_filter_af.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25;DP=30	GT:DP	0/1:30	0/0:25
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
EOF
fi

echo "Test 2: Filtering by AF>=0.2"
$VCFX_EXECUTABLE --filter "AF>=0.2" < data/record_filter_input.vcf > out/record_filter_af.vcf
diff -u expected/record_filter_af.vcf out/record_filter_af.vcf
echo "  Test 2 passed."

# Test 3: Filter by FILTER field
if [ ! -f expected/record_filter_pass.vcf ]; then
  cat > expected/record_filter_pass.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25;DP=30	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	80	PASS	AF=0.10;DP=60	GT:DP	0/0:28	0/1:32
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
EOF
fi

echo "Test 3: Filtering by FILTER==PASS"
$VCFX_EXECUTABLE --filter "FILTER==PASS" < data/record_filter_input.vcf > out/record_filter_pass.vcf
diff -u expected/record_filter_pass.vcf out/record_filter_pass.vcf
echo "  Test 3 passed."

# Test 4: Multiple criteria with AND logic
if [ ! -f expected/record_filter_multiple_and.vcf ]; then
  cat > expected/record_filter_multiple_and.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25;DP=30	GT:DP	0/1:30	0/0:25
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
EOF
fi

echo "Test 4: Multiple criteria with AND logic"
$VCFX_EXECUTABLE --filter "FILTER==PASS;AF>=0.2" --logic and < data/record_filter_input.vcf > out/record_filter_multiple_and.vcf
diff -u expected/record_filter_multiple_and.vcf out/record_filter_multiple_and.vcf
echo "  Test 4 passed."

# Test 5: Multiple criteria with OR logic
if [ ! -f expected/record_filter_multiple_or.vcf ]; then
  cat > expected/record_filter_multiple_or.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	10000	rs123	A	G	100	PASS	AF=0.25;DP=30	GT:DP	0/1:30	0/0:25
chr1	20000	rs456	C	T	80	PASS	AF=0.10;DP=60	GT:DP	0/0:28	0/1:32
chr1	30000	rs789	G	A	50	LowQual	AF=0.05;DP=15	GT:DP	0/1:15	0/0:14
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
EOF
fi

echo "Test 5: Multiple criteria with OR logic"
$VCFX_EXECUTABLE --filter "QUAL>=50;AF>=0.2" --logic or < data/record_filter_input.vcf > out/record_filter_multiple_or.vcf
diff -u expected/record_filter_multiple_or.vcf out/record_filter_multiple_or.vcf
echo "  Test 5 passed."

# Test 6: Memory-mapped file input (positional argument)
echo "Test 6: Memory-mapped file input (positional argument)"
$VCFX_EXECUTABLE --filter "QUAL>=90" data/record_filter_input.vcf > out/record_filter_mmap.vcf
diff -u expected/record_filter_qual.vcf out/record_filter_mmap.vcf
echo "  Test 6 passed."

# Test 7: Memory-mapped file input with -i flag
echo "Test 7: Memory-mapped file input with -i flag"
$VCFX_EXECUTABLE --filter "AF>=0.2" -i data/record_filter_input.vcf > out/record_filter_mmap_i.vcf
diff -u expected/record_filter_af.vcf out/record_filter_mmap_i.vcf
echo "  Test 7 passed."

# Test 8: Multiple criteria with file input
echo "Test 8: Multiple criteria with file input"
$VCFX_EXECUTABLE --filter "FILTER==PASS;AF>=0.2" --logic and data/record_filter_input.vcf > out/record_filter_mmap_multiple.vcf
diff -u expected/record_filter_multiple_and.vcf out/record_filter_mmap_multiple.vcf
echo "  Test 8 passed."

# Test 9: POS filtering
if [ ! -f expected/record_filter_pos.vcf ]; then
  cat > expected/record_filter_pos.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	30000	rs789	G	A	50	LowQual	AF=0.05;DP=15	GT:DP	0/1:15	0/0:14
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
chr1	50000	rs102	G	T	30	LowQual	AF=0.02;DP=10	GT:DP	0/0:10	0/1:8
EOF
fi

echo "Test 9: POS filtering"
$VCFX_EXECUTABLE --filter "POS>=30000" data/record_filter_input.vcf > out/record_filter_pos.vcf
diff -u expected/record_filter_pos.vcf out/record_filter_pos.vcf
echo "  Test 9 passed."

# Test 10: DP (INFO field) filtering
if [ ! -f expected/record_filter_dp.vcf ]; then
  cat > expected/record_filter_dp.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	20000	rs456	C	T	80	PASS	AF=0.10;DP=60	GT:DP	0/0:28	0/1:32
chr1	40000	rs101	T	C	90	PASS	AF=0.35;DP=45	GT:DP	0/1:45	1/1:42
EOF
fi

echo "Test 10: DP (INFO field) filtering"
$VCFX_EXECUTABLE --filter "DP>=40" data/record_filter_input.vcf > out/record_filter_dp.vcf
diff -u expected/record_filter_dp.vcf out/record_filter_dp.vcf
echo "  Test 10 passed."

echo "All VCFX_record_filter tests passed!" 