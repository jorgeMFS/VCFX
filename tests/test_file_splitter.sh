#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
FILE_SPLITTER="$ROOT_DIR/build/src/VCFX_file_splitter/VCFX_file_splitter"

# Create test data directory
TEST_DIR="${SCRIPT_DIR}/data/file_splitter"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Testing VCFX_file_splitter ==="

# Test 1: Basic functionality with multiple chromosomes
echo "Test 1: Basic functionality with multiple chromosomes"
cat > basic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
2	100	rs3	G	A	100	PASS	.
2	200	rs4	T	C	100	PASS	.
X	100	rs5	A	G	100	PASS	.
Y	100	rs6	C	T	100	PASS	.
EOF

# Run the file splitter with default prefix
$FILE_SPLITTER < basic.vcf

# Check if expected output files exist
if [ -f "split_1.vcf" ] && [ -f "split_2.vcf" ] && [ -f "split_X.vcf" ] && [ -f "split_Y.vcf" ]; then
    echo "✓ Output files created successfully"
else
    echo "❌ Test 1 failed: Expected output files not created"
    ls -la
    exit 1
fi

# Check if all files have the same header
for chr in 1 2 X Y; do
    if ! grep -q "##fileformat=VCFv4.2" "split_${chr}.vcf" || \
       ! grep -q "##source=VCFX_test" "split_${chr}.vcf" || \
       ! grep -q "#CHROM" "split_${chr}.vcf"; then
        echo "❌ Test 1 failed: Missing header in split_${chr}.vcf"
        cat "split_${chr}.vcf"
        exit 1
    fi
done

# Check if variants are correctly distributed
if grep -c -v "^#" "split_1.vcf" | grep -q "2" && \
   grep -c -v "^#" "split_2.vcf" | grep -q "2" && \
   grep -c -v "^#" "split_X.vcf" | grep -q "1" && \
   grep -c -v "^#" "split_Y.vcf" | grep -q "1"; then
    echo "✓ Test 1 passed: Variants correctly distributed"
else
    echo "❌ Test 1 failed: Variants not correctly distributed"
    echo "File split_1.vcf:"
    cat "split_1.vcf"
    echo "File split_2.vcf:"
    cat "split_2.vcf"
    echo "File split_X.vcf:"
    cat "split_X.vcf"
    echo "File split_Y.vcf:"
    cat "split_Y.vcf"
    exit 1
fi

# Clean up
rm -f split_*.vcf

# Test 2: Custom prefix
echo "Test 2: Custom prefix"
$FILE_SPLITTER --prefix "chr" < basic.vcf

# Check if expected output files exist with custom prefix
if [ -f "chr_1.vcf" ] && [ -f "chr_2.vcf" ] && [ -f "chr_X.vcf" ] && [ -f "chr_Y.vcf" ]; then
    echo "✓ Test 2 passed: Output files created with custom prefix"
else
    echo "❌ Test 2 failed: Expected output files with custom prefix not created"
    ls -la
    exit 1
fi

# Clean up
rm -f chr_*.vcf

# Test 3: Empty file
echo "Test 3: Empty file"
> empty.vcf

$FILE_SPLITTER < empty.vcf 2> empty.err

if grep -q "Note: No variant data lines were found in the input." empty.err; then
    echo "✓ Test 3 passed: Empty file handled correctly"
else
    echo "❌ Test 3 failed: Empty file not handled correctly"
    cat empty.err
    exit 1
fi

# Test 4: File with only header
echo "Test 4: File with only header"
cat > header_only.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
##contig=<ID=1,length=248956422>
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

$FILE_SPLITTER < header_only.vcf 2> header_only.err

if grep -q "Note: No variant data lines were found in the input." header_only.err; then
    echo "✓ Test 4 passed: Header-only file handled correctly"
else
    echo "❌ Test 4 failed: Header-only file not handled correctly"
    cat header_only.err
    exit 1
fi

# Test 5: Malformed VCF (no CHROM column in data line)
echo "Test 5: Malformed VCF (no CHROM column in data line)"
cat > malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
EOF

$FILE_SPLITTER < malformed.vcf 2> malformed.err

# Check that the file was correctly created for CHROM 100
if [ -f "split_100.vcf" ] && [ -f "split_1.vcf" ]; then
    echo "✓ Test 5 passed: Malformed VCF handled correctly"
else
    echo "❌ Test 5 failed: Malformed VCF not handled correctly"
    ls -la
    exit 1
fi

# Clean up
rm -f split_*.vcf

# Test 6: VCF with additional headers after data lines
echo "Test 6: VCF with additional headers after data lines"
cat > mixed_header.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
##additional_header=value
2	100	rs3	G	A	100	PASS	.
EOF

$FILE_SPLITTER < mixed_header.vcf

# Verify that additional header lines are replicated to all split files
if grep -q "##additional_header=value" "split_1.vcf" && \
   grep -q "##additional_header=value" "split_2.vcf"; then
    echo "✓ Test 6 passed: Additional headers correctly replicated"
else
    echo "❌ Test 6 failed: Additional headers not correctly replicated"
    echo "File split_1.vcf:"
    cat "split_1.vcf"
    echo "File split_2.vcf:"
    cat "split_2.vcf"
    exit 1
fi

# Clean up
rm -f split_*.vcf

# Test 7: Large number of chromosomes
echo "Test 7: Large number of chromosomes"
> large.vcf
echo "##fileformat=VCFv4.2" >> large.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> large.vcf

# Generate variants across 22 chromosomes
for chr in {1..22}; do
    for pos in {1..5}; do
        echo "$chr	$((pos*100))	rs${chr}_${pos}	A	G	100	PASS	." >> large.vcf
    done
done

$FILE_SPLITTER --prefix "large" < large.vcf

# Count the number of created files
file_count=$(ls large_*.vcf | wc -l)
if [ "$file_count" -eq 22 ]; then
    echo "✓ Test 7 passed: Created correct number of files ($file_count)"
else
    echo "❌ Test 7 failed: Expected 22 files, got $file_count"
    ls -la large_*.vcf
    exit 1
fi

# Verify each chromosome file has the expected number of variants
incorrect_files=0
for chr in {1..22}; do
    variant_count=$(grep -v "^#" "large_${chr}.vcf" | wc -l)
    if [ "$variant_count" -ne 5 ]; then
        echo "❌ File large_${chr}.vcf has $variant_count variants instead of 5"
        incorrect_files=$((incorrect_files + 1))
    fi
done

if [ "$incorrect_files" -eq 0 ]; then
    echo "✓ All chromosome files contain the correct number of variants"
else
    echo "❌ Test 7 failed: $incorrect_files files have incorrect variant counts"
    exit 1
fi

# Clean up
rm -f large_*.vcf

# Test 8: Special characters in chromosome names
echo "Test 8: Special characters in chromosome names"
cat > special_chr.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs1	A	G	100	PASS	.
chr2	200	rs2	C	T	100	PASS	.
chr_X	100	rs3	G	A	100	PASS	.
chr.Y	200	rs4	T	C	100	PASS	.
chrM	100	rs5	A	G	100	PASS	.
EOF

$FILE_SPLITTER --prefix "special" < special_chr.vcf

# Check if expected output files exist
if [ -f "special_chr1.vcf" ] && \
   [ -f "special_chr2.vcf" ] && \
   [ -f "special_chr_X.vcf" ] && \
   [ -f "special_chr.Y.vcf" ] && \
   [ -f "special_chrM.vcf" ]; then
    echo "✓ Test 8 passed: Special chromosome names handled correctly"
else
    echo "❌ Test 8 failed: Expected output files for special chromosome names not created"
    ls -la special_*.vcf
    exit 1
fi

# Clean up
rm -f special_*.vcf

# Test 9: VCF with sample columns
echo "Test 9: VCF with sample columns"
cat > samples.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	1/1
2	100	rs3	G	A	100	PASS	.	GT	1/1	0/1
EOF

$FILE_SPLITTER --prefix "samples" < samples.vcf

# Check if sample columns are preserved in output files
if grep -q "SAMPLE1" "samples_1.vcf" && \
   grep -q "SAMPLE2" "samples_1.vcf" && \
   grep -q "GT" "samples_1.vcf" && \
   grep -q "0/0" "samples_1.vcf" && \
   grep -q "0/1" "samples_1.vcf"; then
    echo "✓ Test 9 passed: Sample columns preserved correctly"
else
    echo "❌ Test 9 failed: Sample columns not preserved correctly"
    cat "samples_1.vcf"
    exit 1
fi

# Clean up
rm -f samples_*.vcf

# Test 10: Empty lines in VCF
echo "Test 10: Empty lines in VCF"
cat > empty_lines.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

1	100	rs1	A	G	100	PASS	.

2	100	rs3	G	A	100	PASS	.

EOF

$FILE_SPLITTER --prefix "empty_lines" < empty_lines.vcf

# Check if output files were created successfully despite empty lines
if [ -f "empty_lines_1.vcf" ] && [ -f "empty_lines_2.vcf" ] && \
   grep -q "rs1" "empty_lines_1.vcf" && \
   grep -q "rs3" "empty_lines_2.vcf"; then
    echo "✓ Test 10 passed: Empty lines handled correctly"
else
    echo "❌ Test 10 failed: Empty lines not handled correctly"
    echo "File empty_lines_1.vcf:"
    cat "empty_lines_1.vcf"
    echo "File empty_lines_2.vcf:"
    cat "empty_lines_2.vcf"
    exit 1
fi

# Clean up
rm -f empty_lines_*.vcf

# Test 11: Help message
echo "Test 11: Help message"
$FILE_SPLITTER --help > help_output.txt

if grep -q "VCFX_file_splitter:" help_output.txt && \
   grep -q "Split a VCF file" help_output.txt && \
   grep -q "prefix" help_output.txt; then
    echo "✓ Test 11 passed: Help message is displayed correctly"
else
    echo "❌ Test 11 failed: Help message not displayed correctly"
    cat help_output.txt
    exit 1
fi

echo "All tests for VCFX_file_splitter passed!" 