#!/bin/bash

# Run all tests in the test suite

# Stop on first error
set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# The root directory is one level up from the tests directory
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

# Create build directory if it doesn't exist
mkdir -p "${ROOT_DIR}/build"
cd "${ROOT_DIR}/build"

# Configure and build
cmake .. 
make -j

# Return to the tests directory
cd "${SCRIPT_DIR}"

# List of all test scripts
TEST_SCRIPTS=(
    "test_af_subsetter.sh"
    "test_alignment_checker.sh"
    "test_allele_balance_calc.sh"
    "test_allele_balance_filter.sh"
    "test_allele_counter.sh"
    "test_allele_freq_calc.sh"
    "test_ancestry_assigner.sh"
    "test_ancestry_inferrer.sh"
    "test_annotation_extractor.sh"
    "test_compressor.sh"
    "test_concordance_checker.sh"
    "test_cross_sample_concordance.sh"
    "test_custom_annotator.sh"
    "test_diff_tool.sh"
    "test_distance_calculator.sh"
    "test_dosage_calculator.sh"
    "test_duplicate_remover.sh"
    "test_fasta_converter.sh"
    "test_field_extractor.sh"
    "test_file_splitter.sh"
    "test_format_converter.sh"
    "test_genotype_query.sh"
    "test_gl_filter.sh"
    "test_haplotype_extractor.sh"
    "test_header_parser.sh"
    "test_hwe_tester.sh"
    "test_impact_filter.sh"
    "test_indel_normalizer.sh"
    "test_indexer.sh"
    "test_info_aggregator.sh"
    "test_info_summarizer.sh"
    "test_inbreeding_calculator.sh"
    "test_ld_calculator.sh"
    "test_metadata_summarizer.sh"
    "test_merger.sh"
    "test_missing_data_handler.sh"
    "test_missing_detector.sh"
    "test_multiallelic_splitter.sh"
    "test_nonref_filter.sh"
    "test_outlier_detector.sh"
    "test_phase_checker.sh"
    "test_phase_quality_filter.sh"
    "test_phred_filter.sh"
    "test_population_filter.sh"
    "test_position_subsetter.sh"
    "test_probability_filter.sh"
    "test_quality_adjuster.sh"
    "test_record_filter.sh"
    "test_ref_comparator.sh"
    "test_reformatter.sh"
    "test_region_subsampler.sh"
    "test_sample_extractor.sh"
    "test_sorter.sh"
    "test_sv_handler.sh"
    "test_subsampler.sh"
    "test_validator.sh"
    "test_variant_classifier.sh"
    "test_variant_counter.sh"
    "test_python_bindings.sh"
)

# Run all tests
for TEST_SCRIPT in "${TEST_SCRIPTS[@]}"; do
    echo "---------------------------------------------------------"
    echo "Running $TEST_SCRIPT"
    echo "---------------------------------------------------------"
    
    # Make the script executable if it isn't already
    chmod +x "$TEST_SCRIPT"
    
    # Run the test script
    ./"$TEST_SCRIPT"
    
    echo
done

# Check if Docker is installed and run Docker tests if available
if command -v docker &> /dev/null; then
    echo "---------------------------------------------------------"
    echo "Running Docker tests (test_docker.sh)"
    echo "---------------------------------------------------------"
    
    # Make the Docker test script executable if it isn't already
    chmod +x "test_docker.sh"
    
    # Run the Docker test script
    ./test_docker.sh
    
    echo
else
    echo "---------------------------------------------------------"
    echo "Docker not installed - skipping Docker tests"
    echo "---------------------------------------------------------"
fi

echo "All tests passed!"
