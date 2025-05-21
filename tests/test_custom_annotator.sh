#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

# Create test data directory if it doesn't exist
mkdir -p "$SCRIPT_DIR/data"

# Create test VCF file
cat > "$SCRIPT_DIR/data/test_input.vcf" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	100	PASS	DP=10	GT	0/1
1	200	rs2	T	C,G	100	PASS	DP=20	GT	0/1
1	300	rs3	G	A,T	100	PASS	DP=30	GT	1/2
2	100	rs4	C	T	100	PASS	.	GT	0/1
2	200	rs5	A	G	100	PASS	DP=40	GT	1/1
EOF

# Create annotation file
cat > "$SCRIPT_DIR/data/annotations.txt" << 'EOF'
1	100	A	G	HighImpact
1	200	T	C	ModerateImpact
1	200	T	G	LowImpact
1	300	G	A	NoImpact
2	100	C	T	CriticalImpact
EOF

# Create expected output file
cat > "$SCRIPT_DIR/data/expected_output.vcf" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##INFO=<ID=CustomAnnotation,Number=.,Type=String,Description="Custom annotations added by VCFX_custom_annotator (multi-allelic)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	100	PASS	DP=10;CustomAnnotation=HighImpact	GT	0/1
1	200	rs2	T	C,G	100	PASS	DP=20;CustomAnnotation=ModerateImpact,LowImpact	GT	0/1
1	300	rs3	G	A,T	100	PASS	DP=30;CustomAnnotation=NoImpact,NA	GT	1/2
2	100	rs4	C	T	100	PASS	CustomAnnotation=CriticalImpact	GT	0/1
2	200	rs5	A	G	100	PASS	DP=40;CustomAnnotation=NA	GT	1/1
EOF

echo "=== Testing VCFX_custom_annotator ==="

# Test 1: Basic annotation functionality
echo "Test 1: Basic annotation functionality"
"$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" --add-annotation "$SCRIPT_DIR/data/annotations.txt" < "$SCRIPT_DIR/data/test_input.vcf" > "$SCRIPT_DIR/data/output.vcf"
if diff "$SCRIPT_DIR/data/output.vcf" "$SCRIPT_DIR/data/expected_output.vcf" > /dev/null; then
    echo "✓ Test 1 passed"
else
    echo "✗ Test 1 failed"
    exit 1
fi

# Test 2: Empty annotation file
echo "Test 2: Empty annotation file handling"
touch "$SCRIPT_DIR/data/empty_annotations.txt"
"$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" --add-annotation "$SCRIPT_DIR/data/empty_annotations.txt" < "$SCRIPT_DIR/data/test_input.vcf" > "$SCRIPT_DIR/data/output_empty.vcf"
if grep -q "CustomAnnotation=NA" "$SCRIPT_DIR/data/output_empty.vcf"; then
    echo "✓ Test 2 passed"
else
    echo "✗ Test 2 failed"
    exit 1
fi

# Test 3: Malformed annotation file
echo "Test 3: Malformed annotation file handling"
cat > "$SCRIPT_DIR/data/malformed_annotations.txt" << 'EOF'
1	100	A	G	HighImpact
1	200	T	C	ModerateImpact
1	300	G	A	NoImpact
2	100	C	T	CriticalImpact
malformed line
EOF
"$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" --add-annotation "$SCRIPT_DIR/data/malformed_annotations.txt" < "$SCRIPT_DIR/data/test_input.vcf" > "$SCRIPT_DIR/data/output_malformed.vcf" 2>&1
if grep -q "Warning: Skipping invalid annotation line" "$SCRIPT_DIR/data/output_malformed.vcf"; then
    echo "✓ Test 3 passed"
else
    echo "✗ Test 3 failed"
    exit 1
fi

# Test 4: Missing annotation file
echo "Test 4: Missing annotation file handling"
"$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" --add-annotation "nonexistent.txt" < "$SCRIPT_DIR/data/test_input.vcf" 2>&1 | grep -q "Error: Unable to open annotation file"
if [ $? -eq 0 ]; then
    echo "✓ Test 4 passed"
else
    echo "✗ Test 4 failed"
    exit 1
fi

# Test 5: Help message
echo "Test 5: Help message"
"$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" --help | grep -q "VCFX_custom_annotator: Add custom annotations"
if [ $? -eq 0 ]; then
    echo "✓ Test 5 passed"
else
    echo "✗ Test 5 failed"
    exit 1
fi

# Test 6: No arguments
echo "Test 6: No arguments handling"
"$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" 2>&1 | grep -q "VCFX_custom_annotator: Add custom annotations"
if [ $? -eq 0 ]; then
    echo "✓ Test 6 passed"
else
    echo "✗ Test 6 failed"
    exit 1
fi

# Test 7: Large file handling
echo "Test 7: Large file handling"
# Create a large VCF file
for i in $(seq 1 1000); do
    echo "1	$i	rs$i	A	G	100	PASS	DP=10	GT	0/1"
done > "$SCRIPT_DIR/data/large_input.vcf"
# Create corresponding annotations
for i in $(seq 1 1000); do
    echo "1	$i	A	G	Annotation$i"
done > "$SCRIPT_DIR/data/large_annotations.txt"
# Add VCF header
sed -i '1i\
##fileformat=VCFv4.2\
##contig=<ID=1,length=1000>\
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1\
' "$SCRIPT_DIR/data/large_input.vcf"

time "$ROOT_DIR/build/src/VCFX_custom_annotator/VCFX_custom_annotator" --add-annotation "$SCRIPT_DIR/data/large_annotations.txt" < "$SCRIPT_DIR/data/large_input.vcf" > "$SCRIPT_DIR/data/large_output.vcf"
if [ $? -eq 0 ]; then
    echo "✓ Test 7 passed"
else
    echo "✗ Test 7 failed"
    exit 1
fi

echo "All tests for VCFX_custom_annotator passed!" 

