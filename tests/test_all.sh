#!/usr/bin/env bash
set -e

# Run all tests
echo "Running all VCFX tests..."

# Original tests
bash test_af_subsetter.sh
bash test_alignment_checker.sh
bash test_allele_balance_filter.sh
bash test_allele_balance_calc.sh
bash test_allele_counter.sh

# New tests
bash test_header_parser.sh
bash test_record_filter.sh
bash test_field_extractor.sh
bash test_format_converter.sh
bash test_variant_counter.sh
bash test_sample_extractor.sh

# Additional comprehensive tests
bash test_distance_calculator.sh
bash test_variant_classifier.sh
bash test_cross_sample_concordance.sh
bash test_phase_checker.sh
bash test_validator.sh
bash test_indel_normalizer.sh
bash test_multiallelic_splitter.sh
bash test_sv_handler.sh
bash test_sorter.sh
bash test_position_subsetter.sh

echo "All tests completed successfully!"
