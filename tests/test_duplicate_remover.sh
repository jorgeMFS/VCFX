#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
DUPLICATE_REMOVER="$ROOT_DIR/build/src/VCFX_duplicate_remover/VCFX_duplicate_remover"

# Create test data directory
TEST_DIR="${SCRIPT_DIR}/data/duplicate_remover"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Testing VCFX_duplicate_remover ==="

# Test 1: Basic functionality with exact duplicates
echo "Test 1: Basic functionality with exact duplicates"
cat > exact_duplicates.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
1	100	rs1	A	G	100	PASS	.
2	300	rs3	G	A	100	PASS	.
2	300	rs3	G	A	100	PASS	.
EOF

cat > expected_exact_duplicates.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
2	300	rs3	G	A	100	PASS	.
EOF

$DUPLICATE_REMOVER < exact_duplicates.vcf > output_exact_duplicates.vcf

if diff -q output_exact_duplicates.vcf expected_exact_duplicates.vcf > /dev/null; then
    echo "✓ Test 1 passed"
else
    echo "❌ Test 1 failed"
    echo "Expected:"
    cat expected_exact_duplicates.vcf
    echo "Actual:"
    cat output_exact_duplicates.vcf
    exit 1
fi

# Test 2: Multi-allelic variants with different ordering
echo "Test 2: Multi-allelic variants with different ordering"
cat > multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G,T	100	PASS	.
1	100	rs1	A	T,G	100	PASS	.
1	200	rs2	C	T,G,A	100	PASS	.
1	200	rs2	C	A,G,T	100	PASS	.
1	200	rs2	C	G,T,A	100	PASS	.
EOF

cat > expected_multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G,T	100	PASS	.
1	200	rs2	C	T,G,A	100	PASS	.
EOF

$DUPLICATE_REMOVER < multi_allelic.vcf > output_multi_allelic.vcf

if diff -q output_multi_allelic.vcf expected_multi_allelic.vcf > /dev/null; then
    echo "✓ Test 2 passed"
else
    echo "❌ Test 2 failed"
    echo "Expected:"
    cat expected_multi_allelic.vcf
    echo "Actual:"
    cat output_multi_allelic.vcf
    exit 1
fi

# Test 3: Different entries with same positions but different attributes
echo "Test 3: Different entries with same positions but different attributes"
cat > different_attributes.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	100	rs1	A	T	100	PASS	.
1	100	rs1	C	G	100	PASS	.
1	100	rs1	A	G	90	FAIL	.
EOF

cat > expected_different_attributes.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	100	rs1	A	T	100	PASS	.
1	100	rs1	C	G	100	PASS	.
EOF

$DUPLICATE_REMOVER < different_attributes.vcf > output_different_attributes.vcf

if diff -q output_different_attributes.vcf expected_different_attributes.vcf > /dev/null; then
    echo "✓ Test 3 passed"
else
    echo "❌ Test 3 failed"
    echo "Expected:"
    cat expected_different_attributes.vcf
    echo "Actual:"
    cat output_different_attributes.vcf
    exit 1
fi

# Test 4: No duplicates
echo "Test 4: No duplicates"
cat > no_duplicates.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
2	300	rs3	G	A	100	PASS	.
2	400	rs4	T	C	100	PASS	.
EOF

$DUPLICATE_REMOVER < no_duplicates.vcf > output_no_duplicates.vcf

if diff -q output_no_duplicates.vcf no_duplicates.vcf > /dev/null; then
    echo "✓ Test 4 passed"
else
    echo "❌ Test 4 failed"
    echo "Expected:"
    cat no_duplicates.vcf
    echo "Actual:"
    cat output_no_duplicates.vcf
    exit 1
fi

# Test 5: Malformed lines handling
echo "Test 5: Malformed lines handling"
cat > malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
invalid_line
1	200	rs2	C	T	100	PASS	.
1	invalid	rs3	G	A	100	PASS	.
EOF

# The tool seems to pass through lines with invalid positions as is
cat > expected_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
1	invalid	rs3	G	A	100	PASS	.
EOF

$DUPLICATE_REMOVER < malformed.vcf > output_malformed.vcf 2> malformed.err

if diff -q output_malformed.vcf expected_malformed.vcf > /dev/null; then
    echo "✓ Test 5 passed"
else
    echo "❌ Test 5 failed"
    echo "Expected:"
    cat expected_malformed.vcf
    echo "Actual:"
    cat output_malformed.vcf
    echo "Errors:"
    cat malformed.err
    exit 1
fi

# Test 6: Empty file handling
echo "Test 6: Empty file handling"
> empty.vcf

$DUPLICATE_REMOVER < empty.vcf > output_empty.vcf

if [ -s output_empty.vcf ]; then
    echo "❌ Test 6 failed: Expected empty output file"
    echo "Actual:"
    cat output_empty.vcf
    exit 1
else
    echo "✓ Test 6 passed"
fi

# Test 7: Headers only
echo "Test 7: Headers only"
cat > headers_only.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
##contig=<ID=1,length=248956422>
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

$DUPLICATE_REMOVER < headers_only.vcf > output_headers_only.vcf

if diff -q output_headers_only.vcf headers_only.vcf > /dev/null; then
    echo "✓ Test 7 passed"
else
    echo "❌ Test 7 failed"
    echo "Expected:"
    cat headers_only.vcf
    echo "Actual:"
    cat output_headers_only.vcf
    exit 1
fi

# Test 8: Help message
echo "Test 8: Help message"
$DUPLICATE_REMOVER --help > help_output.txt

if grep -q "VCFX_duplicate_remover" help_output.txt; then
    echo "✓ Test 8 passed"
else
    echo "❌ Test 8 failed"
    echo "Expected help message not found"
    echo "Actual:"
    cat help_output.txt
    exit 1
fi

# Test 9: Variants with samples columns
echo "Test 9: Variants with samples columns"
cat > with_samples.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	1/1
EOF

cat > expected_with_samples.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	1/1
EOF

$DUPLICATE_REMOVER < with_samples.vcf > output_with_samples.vcf

if diff -q output_with_samples.vcf expected_with_samples.vcf > /dev/null; then
    echo "✓ Test 9 passed"
else
    echo "❌ Test 9 failed"
    echo "Expected:"
    cat expected_with_samples.vcf
    echo "Actual:"
    cat output_with_samples.vcf
    exit 1
fi

# Test 10: Stress test with large number of duplicates
echo "Test 10: Stress test with large number of duplicates"
> large.vcf
echo "##fileformat=VCFv4.2" >> large.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> large.vcf

# Generate 1000 variants, half of which are duplicates
for i in {1..500}; do
    echo "1	$i	rs$i	A	G	100	PASS	." >> large.vcf
    # Add a duplicate
    echo "1	$i	rs$i	A	G	100	PASS	." >> large.vcf
    # Add a variant with same position but different attributes
    echo "1	$i	rs${i}_alt	A	T	100	PASS	." >> large.vcf
done

# Count the total number of variants
total_variants=$(grep -v "^#" large.vcf | wc -l | tr -d ' ')

# Run the tool and measure performance
time $DUPLICATE_REMOVER < large.vcf > output_large.vcf

# Count the number of unique variants in the output
unique_variants=$(grep -v "^#" output_large.vcf | wc -l | tr -d ' ')

if [ "$unique_variants" -eq 1000 ]; then
    echo "✓ Test 10 passed: Found 1000 unique variants from $total_variants total"
else
    echo "❌ Test 10 failed: Expected 1000 unique variants, found $unique_variants"
    exit 1
fi

echo "All tests for VCFX_duplicate_remover passed!" 