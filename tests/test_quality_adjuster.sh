#!/usr/bin/env bash

# Test script for VCFX_quality_adjuster
# Tests various quality adjustment transformations including log, sqrt, square, and identity

# Exit on error
set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
BUILD_DIR="../build"
TOOL="${BUILD_DIR}/src/VCFX_quality_adjuster/VCFX_quality_adjuster"
DATA_DIR="data/quality_adjuster"
EXPECTED_DIR="expected/quality_adjuster"
TMP_DIR="tmp"

# Check if build exists and compile if not
if [ ! -e "$TOOL" ]; then
    echo "VCFX_quality_adjuster not found. Building first..."
    cd ..
    mkdir -p build
    cd build
    cmake ..
    make VCFX_quality_adjuster -j
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
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	30	PASS	AF=0.25	GT:DP	0/1:20	0/0:15
1	200	rs2	C	T	0	PASS	AF=0.5	GT:DP	0/1:18	0/1:25
1	300	rs3	G	A	100	PASS	AF=0.1	GT:DP	0/1:30	0/0:20
1	400	rs4	T	C	10	PASS	AF=0.3	GT:DP	0/1:25	1/1:18
1	500	rs5	G	C	.	PASS	AF=0.35	GT:DP	0/0:15	1/1:18
EOF
fi

# Create file with edge cases
if [ ! -e "$DATA_DIR/edge_cases.vcf" ]; then
    echo "Creating edge case VCF file..."
    cat > "$DATA_DIR/edge_cases.vcf" << 'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	0.000001	PASS	.	GT	0/1
1	200	rs2	C	T	1000000	PASS	.	GT	0/1
1	300	rs3	G	A	invalid	PASS	.	GT	0/1
1	400	rs4	T	C		PASS	.	GT	0/1
1	500	rs5	G	C	-10	PASS	.	GT	0/1
EOF
fi

# Create malformed VCF file
if [ ! -e "$DATA_DIR/malformed.vcf" ]; then
    echo "Creating malformed VCF file..."
    cat > "$DATA_DIR/malformed.vcf" << 'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL
1	100	rs1	A	G	30
1	200	rs2	C	T	40
EOF
fi

# Test function
run_test() {
    local test_name=$1
    local transform=$2
    local input_file=$3
    local extra_args=$4
    
    echo "Running test: $test_name"
    
    # Run the tool
    cmd="$TOOL --adjust-qual $transform $extra_args"
    echo "  Command: $cmd"
    
    # Execute command and save output
    if [ -z "$extra_args" ]; then
        cat "$DATA_DIR/$input_file" | $TOOL --adjust-qual "$transform" > "$TMP_DIR/${test_name}_output.vcf" 2> "$TMP_DIR/${test_name}_err.log"
    else
        cat "$DATA_DIR/$input_file" | $TOOL --adjust-qual "$transform" $extra_args > "$TMP_DIR/${test_name}_output.vcf" 2> "$TMP_DIR/${test_name}_err.log"
    fi
    
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

# Test cases
echo "Testing basic transformations..."
run_test "log_transform" "log" "sample.vcf"
run_test "sqrt_transform" "sqrt" "sample.vcf"
run_test "square_transform" "square" "sample.vcf"
run_test "identity_transform" "identity" "sample.vcf"

echo "Testing no-clamp option..."
run_test "log_transform_no_clamp" "log" "sample.vcf" "--no-clamp"
run_test "square_transform_no_clamp" "square" "sample.vcf" "--no-clamp"

echo "Testing edge cases..."
run_test "log_transform_edge" "log" "edge_cases.vcf"
run_test "sqrt_transform_edge" "sqrt" "edge_cases.vcf"
run_test "square_transform_edge" "square" "edge_cases.vcf"

echo "Testing malformed input..."
run_test "malformed_input" "log" "malformed.vcf"

# Test help message
echo "Testing help message..."
"$TOOL" --help > "$TMP_DIR/help_output.txt"
if grep -q "VCFX_quality_adjuster" "$TMP_DIR/help_output.txt" && \
   grep -q "adjust-qual" "$TMP_DIR/help_output.txt"; then
    echo "✅ Test passed: Help message displayed correctly"
else
    echo "❌ Test failed: Help message not displayed correctly"
    echo "Got:"
    cat "$TMP_DIR/help_output.txt"
    exit 1
fi

# Test invalid transformation
echo "Testing invalid transformation..."
if "$TOOL" --adjust-qual "invalid_transform" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/invalid_output.vcf" 2> "$TMP_DIR/invalid_error.log"; then
    echo "❌ Test failed: Tool should exit with error on invalid transformation"
    exit 1
else
    if grep -q "unsupported transformation" "$TMP_DIR/invalid_error.log"; then
        echo "✅ Test passed: Invalid transformation handled correctly"
    else
        echo "❌ Test failed: Invalid transformation error message not displayed correctly"
        echo "Got:"
        cat "$TMP_DIR/invalid_error.log"
        exit 1
    fi
fi

# Test missing required argument
echo "Testing missing required argument..."
"$TOOL" < "$DATA_DIR/sample.vcf" > "$TMP_DIR/missing_arg_output.vcf" 2> "$TMP_DIR/missing_arg_error.log"
# Check if help message is displayed when no arguments are provided
if grep -q "VCFX_quality_adjuster" "$TMP_DIR/missing_arg_output.vcf" && \
   grep -q "adjust-qual" "$TMP_DIR/missing_arg_output.vcf"; then
    echo "✅ Test passed: Help message displayed when no arguments provided"
else
    echo "❌ Test failed: Help message not displayed when no arguments provided"
    echo "Got:"
    cat "$TMP_DIR/missing_arg_output.vcf"
    cat "$TMP_DIR/missing_arg_error.log"
    exit 1
fi

echo "All VCFX_quality_adjuster tests passed! 🎉"
exit 0 