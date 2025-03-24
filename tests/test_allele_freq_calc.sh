#!/usr/bin/env bash

# Test script for VCFX_allele_freq_calc
# Tests various VCF files to verify correct allele frequency calculation

# Exit on error
set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
BUILD_DIR="../build"
TOOL="${BUILD_DIR}/src/VCFX_allele_freq_calc/VCFX_allele_freq_calc"
DATA_DIR="data/allele_freq_calc"
EXPECTED_DIR="expected/allele_freq_calc"
TMP_DIR="tmp"

# Check if build exists and compile if not
if [ ! -e "$TOOL" ]; then
    echo "VCFX_allele_freq_calc not found. Building first..."
    cd ..
    mkdir -p build
    cd build
    cmake ..
    make VCFX_allele_freq_calc -j
    cd "$SCRIPT_DIR"
fi

# Prepare directories
mkdir -p "$DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$TMP_DIR"

# Create test VCF files if they don't exist
if [ ! -e "$DATA_DIR/simple.vcf" ]; then
    echo "Creating simple test VCF file..."
    cat > "$DATA_DIR/simple.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	GT	0/1	0/0	1/1
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0/1	0/1	0/0
1	300	rs3	G	A	50	PASS	AF=0.1	GT	0/0	0/0	1/1
1	400	rs4	T	C	60	PASS	AF=0.75	GT	1/1	1/1	0/1
EOF
fi

if [ ! -e "$DATA_DIR/multiallelic.vcf" ]; then
    echo "Creating multiallelic test VCF file..."
    cat > "$DATA_DIR/multiallelic.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G,T	30	PASS	AF=0.25,0.10	GT	0/1	0/2	1/2
1	200	rs2	C	T,G,A	40	PASS	AF=0.5,0.1,0.2	GT	0/1	2/3	1/2
EOF
fi

if [ ! -e "$DATA_DIR/missing.vcf" ]; then
    echo "Creating test VCF file with missing data..."
    cat > "$DATA_DIR/missing.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	GT	0/1	./.	1/1
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0/.	0/1	./.
EOF
fi

if [ ! -e "$DATA_DIR/phased.vcf" ]; then
    echo "Creating test VCF file with phased genotypes..."
    cat > "$DATA_DIR/phased.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	GT	0|1	0|0	1|1
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0|1	0|1	0|0
EOF
fi

if [ ! -e "$DATA_DIR/complex.vcf" ]; then
    echo "Creating complex test VCF file..."
    cat > "$DATA_DIR/complex.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	GT:AD:DP	0/1:10,5:15	0/0:20,0:20	1/1:0,25:25
1	200	rs2	C	T	40	PASS	AF=0.5	GT:AD:DP	0/1:8,7:15	0/1:10,10:20	0/0:18,2:20
EOF
fi

if [ ! -e "$DATA_DIR/invalid.vcf" ]; then
    echo "Creating invalid test VCF file..."
    cat > "$DATA_DIR/invalid.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
1	100	rs1	A	G	30	PASS	AF=0.25	GT	0/1	0/0	1/1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	200	rs2	C	T	40	PASS	AF=0.5	GT
1	300
EOF
fi

if [ ! -e "$DATA_DIR/no_gt.vcf" ]; then
    echo "Creating test VCF file with no GT field..."
    cat > "$DATA_DIR/no_gt.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	AD:DP	10,5:15	20,0:20	0,25:25
1	200	rs2	C	T	40	PASS	AF=0.5	AD:DP	8,7:15	10,10:20	18,2:20
EOF
fi

# Function to run a test
run_test() {
    local test_name=$1
    local input_file=$2
    
    echo "Running test: $test_name"
    
    # Run the tool
    cat "$DATA_DIR/$input_file" | "$TOOL" > "$TMP_DIR/${test_name}_output.tsv" 2> "$TMP_DIR/${test_name}_err.log"
    
    # Generate expected output for first run
    if [ ! -e "$EXPECTED_DIR/${test_name}.tsv" ]; then
        echo "  Generating expected output for $test_name"
        cp "$TMP_DIR/${test_name}_output.tsv" "$EXPECTED_DIR/${test_name}.tsv"
    fi
    
    # Compare with expected output
    if diff -q "$TMP_DIR/${test_name}_output.tsv" "$EXPECTED_DIR/${test_name}.tsv" > /dev/null; then
        echo "✅ Test passed: $test_name"
    else
        echo "❌ Test failed: $test_name"
        echo "Expected:"
        cat "$EXPECTED_DIR/${test_name}.tsv"
        echo "Got:"
        cat "$TMP_DIR/${test_name}_output.tsv"
        echo "Error log:"
        cat "$TMP_DIR/${test_name}_err.log"
        exit 1
    fi
}

# Function to verify specific frequencies
verify_frequencies() {
    local test_name=$1
    local input_file=$2
    shift 2
    local expected_freqs=("$@")
    
    echo "Running frequency verification test: $test_name"
    
    # Run the tool
    cat "$DATA_DIR/$input_file" | "$TOOL" > "$TMP_DIR/${test_name}_output.tsv" 2> "$TMP_DIR/${test_name}_err.log"
    
    # Skip header line and extract the frequencies
    local actual_freqs=($(tail -n +2 "$TMP_DIR/${test_name}_output.tsv" | cut -f6))
    
    # Check that we have the right number of frequencies
    if [ ${#actual_freqs[@]} -ne ${#expected_freqs[@]} ]; then
        echo "❌ Test failed: $test_name - Expected ${#expected_freqs[@]} variants, got ${#actual_freqs[@]}"
        echo "Output:"
        cat "$TMP_DIR/${test_name}_output.tsv"
        exit 1
    fi
    
    # Check each frequency
    local all_match=true
    for ((i=0; i<${#expected_freqs[@]}; i++)); do
        if [ "${actual_freqs[$i]}" != "${expected_freqs[$i]}" ]; then
            all_match=false
            echo "❌ Frequency mismatch at variant $((i+1)): Expected ${expected_freqs[$i]}, got ${actual_freqs[$i]}"
        fi
    done
    
    if $all_match; then
        echo "✅ Test passed: $test_name - All frequencies match expected values"
    else
        echo "❌ Test failed: $test_name - Frequency mismatches detected"
        echo "Output:"
        cat "$TMP_DIR/${test_name}_output.tsv"
        exit 1
    fi
}

# Test help message
echo "Testing help message..."
"$TOOL" --help > "$TMP_DIR/help_output.txt" 2>&1
if grep -q "VCFX_allele_freq_calc" "$TMP_DIR/help_output.txt" && \
   grep -q "Allele frequency is computed" "$TMP_DIR/help_output.txt"; then
    echo "✅ Test passed: Help message displayed correctly"
else
    echo "❌ Test failed: Help message not displayed correctly"
    echo "Got:"
    cat "$TMP_DIR/help_output.txt"
    exit 1
fi

# Run the tests
run_test "simple" "simple.vcf"
run_test "multiallelic" "multiallelic.vcf"
run_test "missing" "missing.vcf"
run_test "phased" "phased.vcf"
run_test "complex" "complex.vcf"

# Verify known frequencies for simple cases
verify_frequencies "simple_freqs" "simple.vcf" "0.5000" "0.3333" "0.3333" "0.8333"

# Test invalid file cases
echo "Testing with invalid VCF..."
cat "$DATA_DIR/invalid.vcf" | "$TOOL" > "$TMP_DIR/invalid_output.tsv" 2> "$TMP_DIR/invalid_err.log"
if grep -q "Warning:" "$TMP_DIR/invalid_err.log"; then
    echo "✅ Test passed: Tool correctly warned about invalid VCF"
else
    echo "❌ Test failed: Tool should have warned about invalid VCF"
    echo "Error log:"
    cat "$TMP_DIR/invalid_err.log"
    exit 1
fi

# Test file with no GT field
echo "Testing with no GT field..."
cat "$DATA_DIR/no_gt.vcf" | "$TOOL" > "$TMP_DIR/no_gt_output.tsv" 2> "$TMP_DIR/no_gt_err.log"
if grep -q "CHROM" "$TMP_DIR/no_gt_output.tsv" && [ $(wc -l < "$TMP_DIR/no_gt_output.tsv") -eq 1 ]; then
    echo "✅ Test passed: Tool correctly skipped variants with no GT field"
else
    echo "❌ Test failed: Tool should only output header for VCF with no GT field"
    echo "Output:"
    cat "$TMP_DIR/no_gt_output.tsv"
    exit 1
fi

echo "All VCFX_allele_freq_calc tests passed! 🎉"
exit 0 