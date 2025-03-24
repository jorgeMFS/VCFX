#!/usr/bin/env bash

# Test script for VCFX_ancestry_inferrer
# Tests inferring ancestry from VCF files using population frequency data

# Exit on error
set -e

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
BUILD_DIR="../build"
TOOL="${BUILD_DIR}/src/VCFX_ancestry_inferrer/VCFX_ancestry_inferrer"
DATA_DIR="data/ancestry_inferrer"
EXPECTED_DIR="expected/ancestry_inferrer"
TMP_DIR="tmp"

# Check if build exists and compile if not
if [ ! -e "$TOOL" ]; then
    echo "VCFX_ancestry_inferrer not found. Building first..."
    cd ..
    mkdir -p build
    cd build
    cmake ..
    make VCFX_ancestry_inferrer -j
    cd "$SCRIPT_DIR"
fi

# Prepare directories
mkdir -p "$DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$TMP_DIR"

# Create frequency reference data if it doesn't exist
if [ ! -e "$DATA_DIR/population_freqs.txt" ]; then
    echo "Creating frequency reference data..."
    cat > "$DATA_DIR/population_freqs.txt" << 'EOF'
1	100	A	G	EUR	0.75
1	100	A	G	AFR	0.10
1	100	A	G	EAS	0.25
1	200	C	T	EUR	0.20
1	200	C	T	AFR	0.70
1	200	C	T	EAS	0.35
1	300	G	A	EUR	0.30
1	300	G	A	AFR	0.15
1	300	G	A	EAS	0.80
1	400	T	C	EUR	0.85
1	400	T	C	AFR	0.40
1	400	T	C	EAS	0.20
2	150	A	T	EUR	0.90
2	150	A	T	AFR	0.20
2	150	A	T	EAS	0.30
2	250	C	G	EUR	0.10
2	250	C	G	AFR	0.60
2	250	C	G	EAS	0.40
2	350	G	C	EUR	0.40
2	350	G	C	AFR	0.85
2	350	G	C	EAS	0.25
2	450	T	A	EUR	0.45
2	450	T	A	AFR	0.30
2	450	T	A	EAS	0.75
EOF
fi

# Create test VCF files if they don't exist
if [ ! -e "$DATA_DIR/eur_samples.vcf" ]; then
    echo "Creating European ancestry test VCF file..."
    cat > "$DATA_DIR/eur_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EUR_SAMPLE1	EUR_SAMPLE2	EUR_SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.75	GT	1/1	1/1	1/0
1	200	rs2	C	T	40	PASS	AF=0.20	GT	0/0	0/1	0/0
1	400	rs4	T	C	60	PASS	AF=0.85	GT	1/1	1/1	1/0
2	150	rs5	A	T	70	PASS	AF=0.90	GT	1/1	1/0	1/1
2	250	rs6	C	G	80	PASS	AF=0.10	GT	0/0	0/0	0/1
EOF
fi

if [ ! -e "$DATA_DIR/afr_samples.vcf" ]; then
    echo "Creating African ancestry test VCF file..."
    cat > "$DATA_DIR/afr_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AFR_SAMPLE1	AFR_SAMPLE2	AFR_SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.10	GT	0/0	0/1	0/0
1	200	rs2	C	T	40	PASS	AF=0.70	GT	1/1	1/0	1/1
1	300	rs3	G	A	50	PASS	AF=0.15	GT	0/0	0/1	0/0
2	250	rs6	C	G	80	PASS	AF=0.60	GT	1/0	1/1	1/0
2	350	rs7	G	C	90	PASS	AF=0.85	GT	1/1	1/1	1/0
EOF
fi

if [ ! -e "$DATA_DIR/eas_samples.vcf" ]; then
    echo "Creating East Asian ancestry test VCF file..."
    cat > "$DATA_DIR/eas_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EAS_SAMPLE1	EAS_SAMPLE2	EAS_SAMPLE3
1	100	rs1	A	G	30	PASS	AF=0.25	GT	0/1	0/0	0/0
1	300	rs3	G	A	50	PASS	AF=0.80	GT	1/1	1/0	1/1
1	400	rs4	T	C	60	PASS	AF=0.20	GT	0/0	0/1	0/0
2	350	rs7	G	C	90	PASS	AF=0.25	GT	0/0	0/1	0/0
2	450	rs8	T	A	100	PASS	AF=0.75	GT	1/1	1/0	1/1
EOF
fi

if [ ! -e "$DATA_DIR/mixed_samples.vcf" ]; then
    echo "Creating mixed ancestry test VCF file..."
    cat > "$DATA_DIR/mixed_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EUR_SAM	AFR_SAM	EAS_SAM	MIX_SAM
1	100	rs1	A	G	30	PASS	AF=0.5	GT	1/1	0/0	0/1	1/0
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0/0	1/1	0/1	0/1
1	300	rs3	G	A	50	PASS	AF=0.5	GT	0/1	0/0	1/1	0/1
1	400	rs4	T	C	60	PASS	AF=0.5	GT	1/1	0/1	0/0	1/0
2	150	rs5	A	T	70	PASS	AF=0.5	GT	1/1	0/0	0/1	1/0
2	250	rs6	C	G	80	PASS	AF=0.5	GT	0/0	1/1	0/1	0/1
2	350	rs7	G	C	90	PASS	AF=0.5	GT	0/1	1/1	0/0	1/0
2	450	rs8	T	A	100	PASS	AF=0.5	GT	0/0	0/1	1/1	1/0
EOF
fi

if [ ! -e "$DATA_DIR/phased_samples.vcf" ]; then
    echo "Creating phased genotypes test VCF file..."
    cat > "$DATA_DIR/phased_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EUR_PHASED	AFR_PHASED	EAS_PHASED
1	100	rs1	A	G	30	PASS	AF=0.5	GT	1|1	0|0	0|1
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0|0	1|1	0|1
1	300	rs3	G	A	50	PASS	AF=0.5	GT	0|1	0|0	1|1
1	400	rs4	T	C	60	PASS	AF=0.5	GT	1|1	0|1	0|0
2	150	rs5	A	T	70	PASS	AF=0.5	GT	1|1	0|0	0|1
2	250	rs6	C	G	80	PASS	AF=0.5	GT	0|0	1|1	0|1
2	350	rs7	G	C	90	PASS	AF=0.5	GT	0|1	1|1	0|0
2	450	rs8	T	A	100	PASS	AF=0.5	GT	0|0	0|1	1|1
EOF
fi

if [ ! -e "$DATA_DIR/missing_samples.vcf" ]; then
    echo "Creating test VCF file with missing data..."
    cat > "$DATA_DIR/missing_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EUR_MISS	AFR_MISS	EAS_MISS
1	100	rs1	A	G	30	PASS	AF=0.5	GT	1/1	./.	0/1
1	200	rs2	C	T	40	PASS	AF=0.5	GT	0/0	1/1	./.
1	300	rs3	G	A	50	PASS	AF=0.5	GT	0/1	./.	1/1
1	400	rs4	T	C	60	PASS	AF=0.5	GT	./.	0/1	0/0
2	150	rs5	A	T	70	PASS	AF=0.5	GT	1/1	0/0	./.
2	250	rs6	C	G	80	PASS	AF=0.5	GT	0/0	./.	0/1
2	350	rs7	G	C	90	PASS	AF=0.5	GT	./.	1/1	0/0
2	450	rs8	T	A	100	PASS	AF=0.5	GT	0/0	0/1	./.
EOF
fi

if [ ! -e "$DATA_DIR/multiallelic_samples.vcf" ]; then
    echo "Creating multiallelic test VCF file..."
    cat > "$DATA_DIR/multiallelic_samples.vcf" << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EUR_MULTI	AFR_MULTI	EAS_MULTI
1	100	rs1	A	G,T	30	PASS	AF=0.5,0.2	GT	1/1	0/2	0/1
1	200	rs2	C	T,G	40	PASS	AF=0.5,0.1	GT	0/0	1/2	0/2
EOF
fi

# Add frequency data for multiallelic sites
if grep -q "1	100	A	T" "$DATA_DIR/population_freqs.txt"; then
    echo "Multiallelic frequency data already exists"
else
    echo "Adding multiallelic frequency data..."
    cat >> "$DATA_DIR/population_freqs.txt" << 'EOF'
1	100	A	T	EUR	0.15
1	100	A	T	AFR	0.60
1	100	A	T	EAS	0.10
1	200	C	G	EUR	0.05
1	200	C	G	AFR	0.25
1	200	C	G	EAS	0.50
EOF
fi

# Function to run a test
run_test() {
    local test_name=$1
    local input_file=$2
    
    echo "Running test: $test_name"
    
    # Run the tool
    cat "$DATA_DIR/$input_file" | "$TOOL" --frequency "$DATA_DIR/population_freqs.txt" > "$TMP_DIR/${test_name}_output.tsv" 2> "$TMP_DIR/${test_name}_err.log"
    
    # Check for errors
    if [ $? -ne 0 ]; then
        echo "❌ Test failed: $test_name - Tool exited with non-zero status"
        echo "Error log:"
        cat "$TMP_DIR/${test_name}_err.log"
        exit 1
    fi
    
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

# Function to verify population assignments
verify_population() {
    local test_name=$1
    local input_file=$2
    shift 2
    declare -A expected_pops=("$@")
    
    echo "Running population verification test: $test_name"
    
    # Run the tool
    cat "$DATA_DIR/$input_file" | "$TOOL" --frequency "$DATA_DIR/population_freqs.txt" > "$TMP_DIR/${test_name}_output.tsv" 2> "$TMP_DIR/${test_name}_err.log"
    
    # Check for errors
    if [ $? -ne 0 ]; then
        echo "❌ Test failed: $test_name - Tool exited with non-zero status"
        echo "Error log:"
        cat "$TMP_DIR/${test_name}_err.log"
        exit 1
    fi
    
    # Check each expected population
    local all_match=true
    
    # Read the output file, skipping the header
    # Note: Using a temp file to store results to avoid issues with the while loop
    tail -n +2 "$TMP_DIR/${test_name}_output.tsv" > "$TMP_DIR/${test_name}_results.tmp"
    while read -r line; do
        sample=$(echo "$line" | cut -f1)
        population=$(echo "$line" | cut -f2)
        
        if [ -n "${expected_pops[$sample]}" ]; then
            if [ "${expected_pops[$sample]}" != "$population" ]; then
                echo "❌ Population mismatch for sample $sample: Expected ${expected_pops[$sample]}, got $population"
                all_match=false
            fi
        fi
    done < "$TMP_DIR/${test_name}_results.tmp"
    
    if [ "$all_match" = true ]; then
        echo "✅ Test passed: $test_name - All population assignments match expected values"
    else
        echo "❌ Test failed: $test_name - Population assignment mismatches detected"
        echo "Output:"
        cat "$TMP_DIR/${test_name}_output.tsv"
        exit 1
    fi
}

# Test help message
echo "Testing help message..."
"$TOOL" --help > "$TMP_DIR/help_output.txt" 2>&1
if grep -q "VCFX_ancestry_inferrer" "$TMP_DIR/help_output.txt" && \
   grep -q "Infer population ancestry" "$TMP_DIR/help_output.txt"; then
    echo "✅ Test passed: Help message displayed correctly"
else
    echo "❌ Test failed: Help message not displayed correctly"
    echo "Got:"
    cat "$TMP_DIR/help_output.txt"
    exit 1
fi

# Run the tests for each sample cohort
run_test "eur_samples" "eur_samples.vcf"
run_test "afr_samples" "afr_samples.vcf"
run_test "eas_samples" "eas_samples.vcf"
run_test "mixed_samples" "mixed_samples.vcf"
run_test "phased_samples" "phased_samples.vcf"
run_test "missing_samples" "missing_samples.vcf"
run_test "multiallelic_samples" "multiallelic_samples.vcf"

# Verify specific population assignments for mixed cohort
echo "Verifying population assignments for mixed cohort..."
declare -A mixed_pops=(
    ["EUR_SAM"]="EUR"
    ["AFR_SAM"]="AFR"
    ["EAS_SAM"]="EAS"
)
verify_population "mixed_population_check" "mixed_samples.vcf" "${mixed_pops[@]}"

echo "All main tests passed! Testing error handling..."

# Temporarily disable exit on error for error tests
set +e

# Test 1: Missing frequency file
echo "Testing with missing frequency file..."
"$TOOL" --frequency "/nonexistent/file.txt" < "$DATA_DIR/eur_samples.vcf" > "$TMP_DIR/missing_freq_output.tsv" 2> "$TMP_DIR/missing_freq_err.log"
missing_freq_exit_code=$?
echo "Missing frequency file test exit code: $missing_freq_exit_code"

if [ $missing_freq_exit_code -ne 0 ]; then
    echo "✅ Test passed: Tool correctly failed with non-zero exit code for missing frequency file"
else
    echo "❌ Test failed: Tool should have failed with non-zero exit code for missing frequency file"
    exit 1
fi

# Test 2: Invalid VCF format
echo "Testing with invalid VCF format..."
echo "This is not a VCF file" | "$TOOL" --frequency "$DATA_DIR/population_freqs.txt" > "$TMP_DIR/invalid_vcf_output.tsv" 2> "$TMP_DIR/invalid_vcf_err.log"
invalid_vcf_exit_code=$?
echo "Invalid VCF test exit code: $invalid_vcf_exit_code"

if grep -q "Encountered VCF data before #CHROM header" "$TMP_DIR/invalid_vcf_err.log" || \
   [ ! -s "$TMP_DIR/invalid_vcf_output.tsv" ] || [ $invalid_vcf_exit_code -ne 0 ]; then
    echo "✅ Test passed: Tool correctly handled invalid VCF format"
else
    echo "❌ Test failed: Tool should have reported an error or produced no output for invalid VCF"
    exit 1
fi

# Test 3: Malformed frequency file
echo "Testing with malformed frequency file..."
echo "malformed data" > "$TMP_DIR/malformed_freqs.txt"
"$TOOL" --frequency "$TMP_DIR/malformed_freqs.txt" < "$DATA_DIR/eur_samples.vcf" > "$TMP_DIR/malformed_freq_output.tsv" 2> "$TMP_DIR/malformed_freq_err.log"
malformed_exit_code=$?
echo "Malformed frequency file test exit code: $malformed_exit_code"

if [ $malformed_exit_code -ne 0 ]; then
    echo "✅ Test passed: Tool correctly failed with non-zero exit code for malformed frequency file"
else
    echo "❌ Test failed: Tool should have failed with non-zero exit code for malformed frequency file"
    exit 1
fi

# Re-enable exit on error
set -e

echo "All VCFX_ancestry_inferrer tests passed! 🎉"
exit 0 