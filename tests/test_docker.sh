#!/bin/bash
# This script tests the VCFX Docker image using the existing test files from the tests directory
# Docker image to use for the tests. CI may override this when using a locally
# built image.
VCFX_IMAGE="${VCFX_IMAGE:-ghcr.io/jorgemfs/vcfx:latest}"

# Function to check if command succeeded
check_success() {
  if [ $? -ne 0 ]; then
    echo "❌ Error: $1"
    exit 1
  else
    echo "✅ Success: $1"
  fi
}

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
  echo "⚠️ Docker is not installed. Skipping Docker tests."
  exit 0
fi

echo "🧬 Testing VCFX Docker image with official test files..."

# Pull the latest VCFX image
echo "📥 Pulling the latest VCFX Docker image ($VCFX_IMAGE)..."
if docker pull "$VCFX_IMAGE"; then
  check_success "Pulled VCFX Docker image"
else
  echo "⚠️  Unable to pull $VCFX_IMAGE. Building Docker image locally..."
  docker build -t vcfx:local .
  check_success "Built local Docker image"
  VCFX_IMAGE="vcfx:local"
fi

# Get the directory of this script (tests directory)
TESTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "🔍 Using test files from: ${TESTS_DIR}"

# Create temporary output directory in tests/out
TEMP_OUTPUT="${TESTS_DIR}/out/docker_test_output"
mkdir -p "${TEMP_OUTPUT}"
check_success "Created temporary output directory"

# Test 1: List available tools
echo "📋 Listing available VCFX tools..."
docker run --rm $VCFX_IMAGE 'ls -1 /usr/local/bin/VCFX_* | xargs -n1 basename'
check_success "Listed available tools"

# Test 2: Validator test
echo "🔍 Testing VCFX_validator..."
docker run --rm -v "${TESTS_DIR}:/tests" $VCFX_IMAGE 'cat /tests/data/valid.vcf | VCFX_validator'
check_success "Validated valid.vcf file"

# Test 3: Allele frequency calculator test
echo "🧮 Testing VCFX_allele_freq_calc..."
docker run --rm -v "${TESTS_DIR}:/tests" -v "${TEMP_OUTPUT}:/output" \
  $VCFX_IMAGE 'cat /tests/data/allele_freq_calc/test_input.vcf | VCFX_allele_freq_calc > /output/allele_freqs.tsv'
check_success "Calculated allele frequencies"

# Test 4: Sample extractor test
echo "👥 Testing VCFX_sample_extractor..."
docker run --rm -v "${TESTS_DIR}:/tests" -v "${TEMP_OUTPUT}:/output" \
  $VCFX_IMAGE 'cat /tests/data/valid.vcf | VCFX_sample_extractor --samples SAMPLE1 > /output/sample1.vcf'
check_success "Extracted sample"

# Test 5: Variant classifier test
echo "🔬 Testing VCFX_variant_classifier..."
docker run --rm -v "${TESTS_DIR}:/tests" -v "${TEMP_OUTPUT}:/output" \
  $VCFX_IMAGE 'cat /tests/data/classifier_mixed.vcf | VCFX_variant_classifier --append-info > /output/classified.vcf'
check_success "Classified variants"

# Test 6: Testing a pipeline of commands
echo "🔄 Testing a pipeline of VCFX tools..."
docker run --rm -v "${TESTS_DIR}:/tests" -v "${TEMP_OUTPUT}:/output" \
  $VCFX_IMAGE 'cat /tests/data/valid.vcf | VCFX_validator | VCFX_variant_classifier --append-info | VCFX_allele_freq_calc > /output/pipeline_output.tsv'
check_success "Executed pipeline of tools"

echo "🎉 All Docker tests completed successfully!"
echo "📊 Generated files in ${TEMP_OUTPUT}:"
ls -la "${TEMP_OUTPUT}" | sed 's/^/  /'
echo ""
echo "📚 For more information on how to use VCFX with Docker, see the documentation at https://ieeta-pt.github.io/VCFX/docker/"

# Clean up temporary files
echo "🧹 Cleaning up..."
rm -rf "${TEMP_OUTPUT}"
check_success "Cleaned up temporary files" 
