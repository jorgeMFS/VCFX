#!/usr/bin/env bash
set -e

# Create output directories
mkdir -p out data expected

# Set the path to the executable
EXEC="../build/src/VCFX_indel_normalizer/VCFX_indel_normalizer"

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: $EXEC does not exist. Please build the project first."
    exit 1
fi

# Function to run a test
run_test() {
    local test_num=$1
    local description=$2
    local cmd=$3
    local expected=$4
    local output=$5
    local expect_failure=${6:-false}
    
    echo -n "Test $test_num: $description... "
    
    # Run the command
    eval $cmd
    
    # Check if output matches expected
    if diff -q "$expected" "$output" > /dev/null; then
        echo "PASSED"
        return 0
    else
        if [ "$expect_failure" = "true" ]; then
            echo "PASSED (expected failure)"
            return 0
        else
            echo "FAILED"
            echo "Expected:"
            cat "$expected"
            echo "Got:"
            cat "$output"
            echo "Diff:"
            diff "$expected" "$output"
            return 1
        fi
    fi
}

# Create test data

# Basic INDEL test
cat > data/basic_indel.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	ACTG	A	.	PASS	DP=30	GT:AD:DP	0/1:15,15:30
chr1	200	.	A	ATCG	.	PASS	DP=25	GT:AD:DP	0/1:10,15:25
chr1	300	.	AAAAA	A	.	PASS	DP=20	GT:AD:DP	1/1:0,20:20
EOF

# Expected output for basic test
cat > expected/basic_indel_normalized.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	ACTG	A	.	PASS	DP=30	GT:AD:DP	0/1:15,15:30
chr1	200	.	A	ATCG	.	PASS	DP=25	GT:AD:DP	0/1:10,15:25
chr1	300	.	AAAAA	A	.	PASS	DP=20	GT:AD:DP	1/1:0,20:20
EOF

# Multi-allelic test
cat > data/multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	ACTG	A,ACT,AC	.	PASS	DP=40	GT:AD:DP	0/1:10,10,10,10:40
EOF

# Expected output for multi-allelic test
cat > expected/multi_allelic_normalized.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	ACTG	A	.	PASS	DP=40	GT:AD:DP	0/1:10,10,10,10:40
chr1	102	.	TG	T	.	PASS	DP=40	GT:AD:DP	0/1:10,10,10,10:40
chr1	101	.	CTG	C	.	PASS	DP=40	GT:AD:DP	0/1:10,10,10,10:40
EOF

# Left-alignment test (common prefix)
cat > data/left_align.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10	.	GACGTAC	GACATAC	.	PASS	DP=20	GT	0/1
EOF

# Expected output for left-alignment test
cat > expected/left_align_normalized.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	12	.	CGT	CAT	.	PASS	DP=20	GT	0/1
EOF

# Common suffix test
cat > data/common_suffix.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	20	.	ACGTACGT	ATTTACGT	.	PASS	DP=25	GT	0/1
EOF

# Expected output for common suffix test
cat > expected/common_suffix_normalized.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	20	.	ACGT	ATTT	.	PASS	DP=25	GT	0/1
EOF

# Edge cases test (identical REF and ALT after normalization)
cat > data/edge_case.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	30	.	A	A	.	PASS	DP=30	GT	0/1
chr1	40	.	ACGT	ACGT	.	PASS	DP=30	GT	0/1
EOF

# Expected output for edge cases test
cat > expected/edge_case_normalized.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	30	.	A	A	.	PASS	DP=30	GT	0/1
chr1	40	.	ACGT	ACGT	.	PASS	DP=30	GT	0/1
EOF

# Run tests
run_test 1 "Basic INDEL normalization" "$EXEC < data/basic_indel.vcf > out/basic_indel_normalized.vcf" "expected/basic_indel_normalized.vcf" "out/basic_indel_normalized.vcf"

run_test 2 "Multi-allelic variant splitting" "$EXEC < data/multi_allelic.vcf > out/multi_allelic_normalized.vcf" "expected/multi_allelic_normalized.vcf" "out/multi_allelic_normalized.vcf"

run_test 3 "Left alignment (common prefix)" "$EXEC < data/left_align.vcf > out/left_align_normalized.vcf" "expected/left_align_normalized.vcf" "out/left_align_normalized.vcf"

run_test 4 "Common suffix removal" "$EXEC < data/common_suffix.vcf > out/common_suffix_normalized.vcf" "expected/common_suffix_normalized.vcf" "out/common_suffix_normalized.vcf"

run_test 5 "Edge cases (identical REF/ALT)" "$EXEC < data/edge_case.vcf > out/edge_case_normalized.vcf" "expected/edge_case_normalized.vcf" "out/edge_case_normalized.vcf"

run_test 6 "Help message" "$EXEC --help > out/help.txt" "/dev/null" "out/help.txt" "true"

echo "All tests for VCFX_indel_normalizer passed!" 