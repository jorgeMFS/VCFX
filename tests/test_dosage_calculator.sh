#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
DOSAGE_TOOL="$ROOT_DIR/build/src/VCFX_dosage_calculator/VCFX_dosage_calculator"

# Create test data directory
TEST_DIR="${SCRIPT_DIR}/data/dosage_calculator"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Testing VCFX_dosage_calculator ==="

# Test 1: Basic functionality with simple diploid genotypes
echo "Test 1: Basic functionality with simple diploid genotypes"
cat > basic.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	1/1	0/0
EOF

cat > expected_basic.txt << EOF
CHROM	POS	ID	REF	ALT	Dosages
1	100	rs1	A	G	0,1,2
1	200	rs2	C	T	1,2,0
EOF

$DOSAGE_TOOL < basic.vcf > output_basic.txt

if diff -q output_basic.txt expected_basic.txt > /dev/null; then
    echo "✓ Test 1 passed"
else
    echo "❌ Test 1 failed"
    echo "Expected:"
    cat expected_basic.txt
    echo "Actual:"
    cat output_basic.txt
    exit 1
fi

# Test 2: Multi-allelic variants
echo "Test 2: Multi-allelic variants"
cat > multi_allelic.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G,T	100	PASS	.	GT	0/1	0/2	1/2
1	200	rs2	C	T,G,A	100	PASS	.	GT	0/3	1/2	2/3
EOF

cat > expected_multi_allelic.txt << EOF
CHROM	POS	ID	REF	ALT	Dosages
1	100	rs1	A	G,T	1,1,2
1	200	rs2	C	T,G,A	1,2,2
EOF

$DOSAGE_TOOL < multi_allelic.vcf > output_multi_allelic.txt

if diff -q output_multi_allelic.txt expected_multi_allelic.txt > /dev/null; then
    echo "✓ Test 2 passed"
else
    echo "❌ Test 2 failed"
    echo "Expected:"
    cat expected_multi_allelic.txt
    echo "Actual:"
    cat output_multi_allelic.txt
    exit 1
fi

# Test 3: Phased genotypes
echo "Test 3: Phased genotypes"
cat > phased.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0|0	0|1	1|1
1	200	rs2	C	T	100	PASS	.	GT	0|1	1|1	0|0
EOF

cat > expected_phased.txt << EOF
CHROM	POS	ID	REF	ALT	Dosages
1	100	rs1	A	G	0,1,2
1	200	rs2	C	T	1,2,0
EOF

$DOSAGE_TOOL < phased.vcf > output_phased.txt

if diff -q output_phased.txt expected_phased.txt > /dev/null; then
    echo "✓ Test 3 passed"
else
    echo "❌ Test 3 failed"
    echo "Expected:"
    cat expected_phased.txt
    echo "Actual:"
    cat output_phased.txt
    exit 1
fi

# Test 4: Missing genotypes
echo "Test 4: Missing genotypes"
cat > missing.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	./.	0/1	1/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	./.	0/0
1	300	rs3	G	A	100	PASS	.	GT	.	./.	.
EOF

cat > expected_missing.txt << EOF
CHROM	POS	ID	REF	ALT	Dosages
1	100	rs1	A	G	NA,1,2
1	200	rs2	C	T	1,NA,0
1	300	rs3	G	A	NA,NA,NA
EOF

$DOSAGE_TOOL < missing.vcf > output_missing.txt

if diff -q output_missing.txt expected_missing.txt > /dev/null; then
    echo "✓ Test 4 passed"
else
    echo "❌ Test 4 failed"
    echo "Expected:"
    cat expected_missing.txt
    echo "Actual:"
    cat output_missing.txt
    exit 1
fi

# Test 5: Malformed genotypes
echo "Test 5: Malformed genotypes"
cat > malformed.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/X	0/1	1/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	ABC	0/0
1	300	rs3	G	A	100	PASS	.	GT	0/0/1	0/1	1/1
EOF

cat > expected_malformed.txt << EOF
CHROM	POS	ID	REF	ALT	Dosages
1	100	rs1	A	G	NA,1,2
1	200	rs2	C	T	1,NA,0
1	300	rs3	G	A	NA,1,2
EOF

$DOSAGE_TOOL < malformed.vcf > output_malformed.txt

if diff -q output_malformed.txt expected_malformed.txt > /dev/null; then
    echo "✓ Test 5 passed"
else
    echo "❌ Test 5 failed"
    echo "Expected:"
    cat expected_malformed.txt
    echo "Actual:"
    cat output_malformed.txt
    exit 1
fi

# Test 6: Missing GT field in FORMAT
echo "Test 6: Missing GT field in FORMAT"
cat > missing_gt.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	DP:GQ	10:30	15:40	20:50
1	200	rs2	C	T	100	PASS	.	GT:DP	0/1:10	1/1:15	0/0:20
EOF

cat > expected_missing_gt.txt << EOF
CHROM	POS	ID	REF	ALT	Dosages
1	100	rs1	A	G	NA
1	200	rs2	C	T	1,2,0
EOF

$DOSAGE_TOOL < missing_gt.vcf > output_missing_gt.txt

if diff -q output_missing_gt.txt expected_missing_gt.txt > /dev/null; then
    echo "✓ Test 6 passed"
else
    echo "❌ Test 6 failed"
    echo "Expected:"
    cat expected_missing_gt.txt
    echo "Actual:"
    cat output_missing_gt.txt
    exit 1
fi

# Test 7: Missing CHROM header
echo "Test 7: Missing CHROM header"
cat > missing_header.vcf << EOF
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	1/1	0/0
EOF

# Temporarily disable exit on error for this test
set +e
$DOSAGE_TOOL < missing_header.vcf > output_missing_header.txt 2> error_missing_header.txt
exit_code=$?
set -e

# We expect the tool to report an error about missing header
if grep -q "Error: VCF header (#CHROM) not found before variant records" error_missing_header.txt; then
    echo "✓ Test 7 passed"
else
    echo "❌ Test 7 failed"
    echo "Expected error about missing header"
    echo "Actual error:"
    cat error_missing_header.txt
    exit 1
fi

# Test 8: Help message
echo "Test 8: Help message"
$DOSAGE_TOOL --help > help_output.txt

if grep -q "VCFX_dosage_calculator: Calculate genotype dosage for each variant in a VCF file" help_output.txt; then
    echo "✓ Test 8 passed"
else
    echo "❌ Test 8 failed"
    echo "Expected help message not found"
    echo "Actual output:"
    cat help_output.txt
    exit 1
fi

# Test 9: Performance with large file
echo "Test 9: Performance with large file"

# Create a large VCF file with many variants and samples
> large.vcf
echo "##fileformat=VCFv4.2" >> large.vcf
echo -n "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" >> large.vcf
for i in {1..20}; do
    echo -n "	SAMPLE$i" >> large.vcf
done
echo "" >> large.vcf

for i in {1..1000}; do
    genotypes=""
    for j in {1..20}; do
        # Generate random genotypes (0/0, 0/1, 1/1)
        gt=$((RANDOM % 3))
        case $gt in
            0) genotypes="$genotypes	0/0" ;;
            1) genotypes="$genotypes	0/1" ;;
            2) genotypes="$genotypes	1/1" ;;
        esac
    done
    echo "1	$i	rs$i	A	G	100	PASS	.	GT$genotypes" >> large.vcf
done

# Time the execution
time $DOSAGE_TOOL < large.vcf > output_large.txt

# Check that output has the right number of lines (header + 1000 variants)
line_count=$(wc -l < output_large.txt)
if [ "$line_count" -eq 1001 ]; then
    echo "✓ Test 9 passed"
else
    echo "❌ Test 9 failed"
    echo "Expected 1001 lines (header + 1000 variants)"
    echo "Actual line count: $line_count"
    exit 1
fi

echo "All tests for VCFX_dosage_calculator passed!" 