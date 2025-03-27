#!/usr/bin/env bash

# Test script for VCFX_probability_filter
# Tests various probability filter conditions with different operators

# Exit on error
set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
BUILD_DIR="../build"
TOOL="${BUILD_DIR}/src/VCFX_probability_filter/VCFX_probability_filter"
DATA_DIR="data/probability_filter"
EXPECTED_DIR="expected/probability_filter"
TMP_DIR="tmp"

# Check if build exists and compile if not
if [ ! -e "$TOOL" ]; then
    echo "VCFX_probability_filter not found. Building first..."
    cd ..
    mkdir -p build
    cd build
    cmake ..
    make VCFX_probability_filter -j
    cd "$SCRIPT_DIR"
fi

# Prepare directories
mkdir -p "$DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$TMP_DIR"

# Create test VCF file if it doesn't exist
if [ ! -e "$DATA_DIR/sample.vcf" ]; then
    echo "Creating test VCF file..."
    cat > "$DATA_DIR/sample.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype Probabilities">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	GT:GP	0/1:0.01,0.98,0.01	0/0:0.99,0.01,0	1/1:0,0.02,0.98
1	200	rs2	C	T	40	PASS	AF=0.5	GT:GP	0/1:0.05,0.9,0.05	0/1:0.1,0.8,0.1	0/0:0.95,0.04,0.01
1	300	rs3	G	A	50	PASS	AF=0.1	GT:GP	0/0:0.85,0.15,0	0/0:0.92,0.08,0	1/1:0,0.05,0.95
1	400	rs4	T	C	60	PASS	AF=0.3	GT:GP	0/1:0.1,0.7,0.2	1/1:0,0.1,0.9	0/1:0.1,0.7,0.2
1	500	rs5	G	C	70	PASS	AF=0.35	GT:GP	0/0:0.94,0.05,0.01	1/1:0.01,0.05,0.94	0/1:0.2,0.75,0.05
EOF
fi

# Create missing field test file
if [ ! -e "$DATA_DIR/missing_field.vcf" ]; then
    echo "Creating missing field test file..."
    cat > "$DATA_DIR/missing_field.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	30	PASS	AF=0.25	GT	0/1
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0/0
EOF
fi

# Create missing value test file
if [ ! -e "$DATA_DIR/missing_value.vcf" ]; then
    echo "Creating missing value test file..."
    cat > "$DATA_DIR/missing_value.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype Probabilities">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	30	PASS	AF=0.25	GT:GP	0/1:0.01,0.98,0.01	0/0:0.99,0.01,0
1	200	rs2	C	T	40	PASS	AF=0.5	GT:GP	0/1:0.05,0.9,0.05	0/1:.,.,.
1	300	rs3	G	A	50	PASS	AF=0.1	GT:GP	0/0:.	0/0:0.92,0.08,0
EOF
fi

# Create invalid data file
if [ ! -e "$DATA_DIR/invalid.vcf" ]; then
    echo "Creating invalid data file..."
    cat > "$DATA_DIR/invalid.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype Probabilities">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	30	PASS	AF=0.25	GT:GP	0/1:invalid
1	200	rs2	C	T	40	PASS	AF=0.5	GT:GP	0/1:0.05,0.9,xxx
EOF
fi

# Test function
run_test() {
    local test_name=$1
    local condition=$2
    local input_file=$3
    
    echo "Running test: $test_name"
    
    # Run the tool
    cmd="$TOOL --filter-probability \"$condition\""
    echo "  Command: $cmd"
    
    # Execute command and save output
    cat "$DATA_DIR/$input_file" | $TOOL --filter-probability "$condition" > "$TMP_DIR/${test_name}_output.vcf" 2> "$TMP_DIR/${test_name}_err.log"
    
    # Generate expected output for first run
    if [ ! -e "$EXPECTED_DIR/${test_name}.vcf" ]; then
        echo "  Generating expected output for $test_name"
        cp "$TMP_DIR/${test_name}_output.vcf" "$EXPECTED_DIR/${test_name}.vcf"
    fi
    
    # Compare with expected output
    if diff -q "$TMP_DIR/${test_name}_output.vcf" "$EXPECTED_DIR/${test_name}.vcf" > /dev/null; then
        echo "✅ Test passed: $test_name"
    else
        echo "❌ Test failed: $test_name"
        echo "Expected:"
        cat "$EXPECTED_DIR/${test_name}.vcf"
        echo "Got:"
        cat "$TMP_DIR/${test_name}_output.vcf"
        exit 1
    fi
}

# Test case for error handling
test_error_handling() {
    local test_name=$1
    local condition=$2
    local input_file=$3
    local expected_error=$4
    
    echo "Running error test: $test_name"
    
    # Run the tool
    cmd="$TOOL --filter-probability \"$condition\""
    echo "  Command: $cmd"
    
    # Execute command and save output
    cat "$DATA_DIR/$input_file" | $TOOL --filter-probability "$condition" > "$TMP_DIR/${test_name}_output.vcf" 2> "$TMP_DIR/${test_name}_err.log"
    
    # Check for expected error message
    if grep -q "$expected_error" "$TMP_DIR/${test_name}_err.log"; then
        echo "✅ Test passed: $test_name - Found expected error message"
    else
        echo "❌ Test failed: $test_name - Expected error message not found"
        echo "Expected error containing: $expected_error"
        echo "Got:"
        cat "$TMP_DIR/${test_name}_err.log"
        exit 1
    fi
}

# Test cases for various operators
echo "Testing basic operators..."
run_test "gt_operator" "GP>0.9" "sample.vcf"
run_test "lt_operator" "GP<0.1" "sample.vcf"
run_test "ge_operator" "GP>=0.9" "sample.vcf"
run_test "le_operator" "GP<=0.1" "sample.vcf"
run_test "eq_operator" "GP==0.9" "sample.vcf"
run_test "ne_operator" "GP!=0.9" "sample.vcf"

# Test error handling
echo "Testing error handling..."
test_error_handling "missing_field" "GP>0.9" "missing_field.vcf" "Specified field \"GP\" not found in FORMAT column"
test_error_handling "invalid_condition" "GPX0.9" "sample.vcf" "Invalid filter condition format"

# Test help message
echo "Testing help message..."
"$TOOL" --help > "$TMP_DIR/help_output.txt" 2>&1
if grep -q "VCFX_probability_filter" "$TMP_DIR/help_output.txt" && \
   grep -q "filter-probability" "$TMP_DIR/help_output.txt"; then
    echo "✅ Test passed: Help message displayed correctly"
else
    echo "❌ Test failed: Help message not displayed correctly"
    echo "Got:"
    cat "$TMP_DIR/help_output.txt"
    exit 1
fi

echo "All VCFX_probability_filter tests passed! 🎉"
exit 0 