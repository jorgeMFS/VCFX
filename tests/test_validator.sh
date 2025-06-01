#!/usr/bin/env bash

# Create output directories
mkdir -p out data expected

# Set the path to the executable
EXEC="../build/src/VCFX_validator/VCFX_validator"

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: $EXEC does not exist. Please build the project first."
    exit 1
fi

# Function to run a test with expected success
run_test_success() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local opts="$4"
    
    echo -n "Test $test_num: $description... "
    
    # Run the command using process substitution
    local output
    local exit_code
    output=$($EXEC $opts < "$input_file" 2>&1)
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "PASSED"
        return 0
    else
        echo "FAILED (Expected success, got failure)"
        echo "Exit code: $exit_code"
        echo "Command output: $output"
        return 1
    fi
}

# Function to run a test with expected failure
run_test_failure() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local expected_error=$4
    local opts="$5"
    
    echo -n "Test $test_num: $description... "
    
    # Run the command using process substitution
    local output
    local exit_code
    output=$($EXEC $opts < "$input_file" 2>&1)
    exit_code=$?
    
    if [ $exit_code -ne 0 ]; then
        # Check if the output contains the expected error message
        if echo "$output" | grep -q "$expected_error"; then
            echo "PASSED (Expected failure)"
            return 0
        else
            echo "FAILED (Expected error message not found)"
            echo "Expected to contain: $expected_error"
            echo "Exit code: $exit_code"
            echo "Actual output: $output"
            return 1
        fi
    else
        echo "FAILED (Expected failure, got success)"
        echo "Exit code: $exit_code"
        echo "Command output: $output"
        return 1
    fi
}

# Create test files (same as in test_validator.sh)
# Valid VCF file - basic
cat > data/valid.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10	GT:GQ:DP	0/1:48:8	0/0:43:5
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15	GT:GQ:DP	0/1:56:12	1/1:60:9
chr3	300	.	T	A	90	PASS	NS=2;DP=20	GT:GQ:DP	0/0:50:10	0/1:65:10
EOF

# Invalid VCF - missing CHROM header
cat > data/no_chrom_header.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15
EOF

# Invalid VCF - malformed meta line
cat > data/malformed_meta.vcf << EOF
##fileformat=VCFv4.2
#INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10
EOF

# Invalid VCF - data before CHROM header
cat > data/data_before_header.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15
EOF

# Invalid VCF - too few columns
cat > data/too_few_columns.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
chr1	100	rs123	A	T	60	PASS
EOF

# Invalid VCF - invalid POS (negative)
cat > data/invalid_pos.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	-10	rs123	A	T	60	PASS	NS=2;DP=10
EOF

# Invalid VCF - empty REF
cat > data/empty_ref.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123		T	60	PASS	NS=2;DP=10
EOF

# Invalid VCF - empty ALT
cat > data/empty_alt.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF		QUAL	FILTER	INFO
chr1	100	rs123	A		60	PASS	NS=2;DP=10
EOF

# Invalid VCF - negative QUAL
cat > data/negative_qual.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	T	-60	PASS	NS=2;DP=10
EOF

# Invalid VCF - empty FILTER
cat > data/empty_filter.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL		INFO
chr1	100	rs123	A	T	60		NS=2;DP=10
EOF

# Valid VCF - with missing values represented by dots
cat > data/valid_with_dots.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	.	PASS	.
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15
EOF

# Header has one sample column but a data line includes two sample columns
cat > data/mismatched_columns.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/mismatched_columns.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/0\n' >> data/mismatched_columns.vcf

# FORMAT expects two entries but sample has three
cat > data/format_mismatch.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/format_mismatch.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:DP\t0/1:30:7\n' >> data/format_mismatch.vcf

# Undefined INFO field
cat > data/undefined_info.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> data/undefined_info.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\tDP=5;FOO=1\n' >> data/undefined_info.vcf

# Undefined FORMAT field
cat > data/undefined_format.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/undefined_format.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:XY\t0/1:10\n' >> data/undefined_format.vcf

# Invalid genotype string
cat > data/invalid_genotype.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/invalid_genotype.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0//1\n' >> data/invalid_genotype.vcf

# Invalid REF/ALT bases
cat > data/invalid_bases.vcf << EOF
##fileformat=VCFv4.2
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> data/invalid_bases.vcf
printf 'chr1\t100\t.\tA\tR\t60\tPASS\t.\n' >> data/invalid_bases.vcf

# Duplicate records
cat > data/duplicate_records.vcf << EOF
##fileformat=VCFv4.2
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> data/duplicate_records.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\n' >> data/duplicate_records.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\n' >> data/duplicate_records.vcf

# gzipped input
if [ ! -f data/valid.vcf.gz ]; then
    gzip -c data/valid.vcf > data/valid.vcf.gz
fi

# Run each test separately and track failures
failures=0

# Test 1
run_test_success 1 "Basic valid VCF" "data/valid.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 2
run_test_failure 2 "Missing CHROM header" "data/no_chrom_header.vcf" "data line encountered before #CHROM"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 3
run_test_failure 3 "Malformed meta line" "data/malformed_meta.vcf" "is a header line but neither starts with '##' nor is a #CHROM"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 4
run_test_failure 4 "Data before CHROM header" "data/data_before_header.vcf" "data line encountered before #CHROM"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 5
run_test_failure 5 "Too few columns" "data/too_few_columns.vcf" "has <8 columns"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 6
run_test_failure 6 "Invalid position" "data/invalid_pos.vcf" "POS must be >0"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 7
run_test_failure 7 "Empty REF" "data/empty_ref.vcf" "REF is empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 8
run_test_failure 8 "Empty ALT" "data/empty_alt.vcf" "ALT is empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 9
run_test_failure 9 "Negative QUAL" "data/negative_qual.vcf" "negative QUAL"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 10
run_test_failure 10 "Empty FILTER" "data/empty_filter.vcf" "FILTER is empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 11
run_test_success 11 "Valid VCF with dots" "data/valid_with_dots.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 12 - Help message
echo -n "Test 12: Help message... "
$EXEC --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
    failures=$((failures + 1))
fi

# Test 13 - strict mode valid file
run_test_success 13 "Strict valid VCF" "data/valid.vcf" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 14 - mismatched columns in strict mode
run_test_failure 14 "Strict mismatched columns" "data/mismatched_columns.vcf" "columns" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 15 - FORMAT/sample mismatch in strict mode
run_test_failure 15 "Strict format mismatch" "data/format_mismatch.vcf" "FORMAT" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 16 - undefined INFO field
run_test_failure 16 "Undefined INFO" "data/undefined_info.vcf" "INFO field" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 17 - undefined FORMAT field
run_test_failure 17 "Undefined FORMAT" "data/undefined_format.vcf" "FORMAT field" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 18 - invalid genotype
run_test_failure 18 "Invalid genotype" "data/invalid_genotype.vcf" "invalid genotype" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 19 - invalid bases
run_test_failure 19 "Invalid bases" "data/invalid_bases.vcf" "invalid characters" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 20 - duplicate detection
run_test_failure 20 "Duplicate records" "data/duplicate_records.vcf" "duplicate variant" "--strict --report-dups"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 21 - gzipped input
run_test_success 21 "Gzip input" "data/valid.vcf.gz"
[ $? -ne 0 ] && failures=$((failures + 1))

if [ $failures -eq 0 ]; then
    echo "All tests for VCFX_validator passed!"
    exit 0
else
    echo "$failures test(s) failed."
    exit 1
fi 