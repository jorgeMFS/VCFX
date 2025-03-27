#!/bin/bash

# Test script for VCFX_annotation_extractor

# Stop on first error
set -e

# Change to the directory containing this script
cd "$(dirname "$0")"

# Get the root directory
ROOT_DIR="$(cd .. && pwd)"

# Make sure the build directory exists and the tool is built
if [ ! -f "${ROOT_DIR}/build/src/VCFX_annotation_extractor/VCFX_annotation_extractor" ]; then
    echo "Building VCFX_annotation_extractor..."
    mkdir -p "${ROOT_DIR}/build"
    cd "${ROOT_DIR}/build"
    cmake ..
    make VCFX_annotation_extractor
    cd "${ROOT_DIR}/tests"
fi

# Create a temporary directory for test outputs
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Create test data directory if it doesn't exist
mkdir -p data

# Create test VCF file with various annotation scenarios
cat > data/annotation_input.vcf << 'EOF'
##fileformat=VCFv4.2
##source=VCFX_test
##INFO=<ID=ANN,Number=A,Type=String,Description="Annotation">
##INFO=<ID=Gene,Number=1,Type=String,Description="Gene name">
##INFO=<ID=Impact,Number=1,Type=String,Description="Impact level">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs1	A	G	100	PASS	ANN=missense_variant;Gene=BRCA1;Impact=HIGH;DP=100
chr1	200	rs2	C	T,G	80	PASS	ANN=stop_gained,synonymous_variant;Gene=TP53;Impact=HIGH;DP=50
chr1	300	rs3	G	A,T,C	60	PASS	ANN=missense_variant,stop_gained,intergenic_region;Gene=EGFR;Impact=MODERATE;DP=75
chr2	400	rs4	T	C	90	PASS	Gene=BRCA2;Impact=LOW;DP=200
chr2	500	rs5	A	G,T	70	PASS	ANN=missense_variant,stop_gained;Gene=PTEN;Impact=HIGH;DP=150
EOF

# Create expected output for basic annotation extraction
cat > expected/annotation_extractor_basic.txt << 'EOF'
CHROM	POS	ID	REF	ALT	ANN	Gene	Impact	DP
chr1	100	rs1	A	G	missense_variant	BRCA1	HIGH	100
chr1	200	rs2	C	T	stop_gained	TP53	HIGH	50
chr1	200	rs2	C	G	synonymous_variant	TP53	HIGH	50
chr1	300	rs3	G	A	missense_variant	EGFR	MODERATE	75
chr1	300	rs3	G	T	stop_gained	EGFR	MODERATE	75
chr1	300	rs3	G	C	intergenic_region	EGFR	MODERATE	75
chr2	400	rs4	T	C	NA	BRCA2	LOW	200
chr2	500	rs5	A	G	missense_variant	PTEN	HIGH	150
chr2	500	rs5	A	T	stop_gained	PTEN	HIGH	150
EOF

# Test 1: Basic annotation extraction
echo "Test 1: Basic annotation extraction"
${ROOT_DIR}/build/src/VCFX_annotation_extractor/VCFX_annotation_extractor --annotation-extract "ANN,Gene,Impact,DP" < data/annotation_input.vcf > ${TEMP_DIR}/output.txt

# Compare with expected output
if diff -w -q ${TEMP_DIR}/output.txt expected/annotation_extractor_basic.txt > /dev/null; then
    echo "  ✓ Basic annotation extraction test passed"
else
    echo "  ✗ Basic annotation extraction test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/output.txt expected/annotation_extractor_basic.txt
    exit 1
fi

# Test 2: Single annotation extraction
echo "Test 2: Single annotation extraction"
cat > expected/annotation_extractor_single.txt << 'EOF'
CHROM	POS	ID	REF	ALT	Gene
chr1	100	rs1	A	G	BRCA1
chr1	200	rs2	C	T	TP53
chr1	200	rs2	C	G	TP53
chr1	300	rs3	G	A	EGFR
chr1	300	rs3	G	T	EGFR
chr1	300	rs3	G	C	EGFR
chr2	400	rs4	T	C	BRCA2
chr2	500	rs5	A	G	PTEN
chr2	500	rs5	A	T	PTEN
EOF

${ROOT_DIR}/build/src/VCFX_annotation_extractor/VCFX_annotation_extractor --annotation-extract "Gene" < data/annotation_input.vcf > ${TEMP_DIR}/output.txt

if diff -w -q ${TEMP_DIR}/output.txt expected/annotation_extractor_single.txt > /dev/null; then
    echo "  ✓ Single annotation extraction test passed"
else
    echo "  ✗ Single annotation extraction test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/output.txt expected/annotation_extractor_single.txt
    exit 1
fi

# Test 3: Missing annotation handling
echo "Test 3: Missing annotation handling"
cat > expected/annotation_extractor_missing.txt << 'EOF'
CHROM	POS	ID	REF	ALT	MissingAnn
chr1	100	rs1	A	G	NA
chr1	200	rs2	C	T	NA
chr1	200	rs2	C	G	NA
chr1	300	rs3	G	A	NA
chr1	300	rs3	G	T	NA
chr1	300	rs3	G	C	NA
chr2	400	rs4	T	C	NA
chr2	500	rs5	A	G	NA
chr2	500	rs5	A	T	NA
EOF

${ROOT_DIR}/build/src/VCFX_annotation_extractor/VCFX_annotation_extractor --annotation-extract "MissingAnn" < data/annotation_input.vcf > ${TEMP_DIR}/output.txt

if diff -w -q ${TEMP_DIR}/output.txt expected/annotation_extractor_missing.txt > /dev/null; then
    echo "  ✓ Missing annotation handling test passed"
else
    echo "  ✗ Missing annotation handling test failed"
    echo "Differences:"
    diff -w ${TEMP_DIR}/output.txt expected/annotation_extractor_missing.txt
    exit 1
fi

# Test 4: Error handling - no annotations specified
echo "Test 4: Error handling - no annotations specified"
set +e  # Don't exit on error for this test
${ROOT_DIR}/build/src/VCFX_annotation_extractor/VCFX_annotation_extractor > ${TEMP_DIR}/error_output.txt 2>&1
EXIT_CODE=$?
set -e  # Restore exit on error

if [ $EXIT_CODE -ne 0 ] && grep -q "VCFX_annotation_extractor: Extract variant annotations" ${TEMP_DIR}/error_output.txt; then
    echo "  ✓ Error handling test passed"
else
    echo "  ✗ Error handling test failed"
    echo "Expected non-zero exit code and help message"
    echo "Got exit code $EXIT_CODE and output:"
    cat ${TEMP_DIR}/error_output.txt
    exit 1
fi

# Test 5: Help message
echo "Test 5: Help message"
set +e  # Don't exit on error for help message test
${ROOT_DIR}/build/src/VCFX_annotation_extractor/VCFX_annotation_extractor --help > ${TEMP_DIR}/help_output.txt
set -e  # Restore exit on error

if grep -q "VCFX_annotation_extractor: Extract variant annotations from a VCF file" ${TEMP_DIR}/help_output.txt; then
    echo "  ✓ Help message test passed"
else
    echo "  ✗ Help message test failed"
    echo "Expected help message containing tool description"
    echo "Got output:"
    cat ${TEMP_DIR}/help_output.txt
    exit 1
fi

echo "All tests for VCFX_annotation_extractor passed!" 