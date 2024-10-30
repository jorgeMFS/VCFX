#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define an array of tool names to compile
tools=(
    VCFX_header_parser
    VCFX_record_filter
    VCFX_field_extractor
    VCFX_format_converter
    VCFX_variant_counter
    VCFX_sample_extractor
    VCFX_sorter
    VCFX_validator
    VCFX_subsampler
    VCFX_genotype_query
    VCFX_allele_freq_calc
    VCFX_indexer
    VCFX_compressor
    VCFX_position_subsetter
    VCFX_haplotype_extractor
    VCFX_info_parser
    VCFX_variant_classifier
    VCFX_duplicate_remover
    VCFX_info_summarizer
    VCFX_distance_calculator
    VCFX_multiallelic_splitter
    VCFX_missing_data_handler
    VCFX_concordance_checker
    VCFX_allele_balance_calc
    VCFX_allele_counter
    VCFX_phase_checker
    VCFX_annotation_extractor
    VCFX_phred_filter
    VCFX_merger
    VCFX_metadata_summarizer
    VCFX_hwe_tester
    VCFX_fasta_converter
    VCFX_nonref_filter
    VCFX_dosage_calculator
    VCFX_population_filter
    VCFX_file_splitter
    VCFX_gl_filter
    VCFX_ref_comparator
    VCFX_ancestry_inferrer
    VCFX_impact_filter
    VCFX_info_aggregator
    VCFX_probability_filter
    VCFX_diff_tool
    VCFX_cross_sample_concordance
    VCFX_phase_quality_filter
    VCFX_indel_normalizer
    VCFX_custom_annotator
    VCFX_region_subsampler
    VCFX_allele_balance_filter
    VCFX_missing_detector
    VCFX_haplotype_phaser
    VCFX_af_subsetter
    VCFX_sv_handler
    VCFX_reformatter
    VCFX_quality_adjuster
    VCFX_inbreeding_calculator
    VCFX_outlier_detector
    VCFX_alignment_checker
    VCFX_ancestry_assigner
    VCFX_ld_calculator
)

# Output directory for WebAssembly files
output_dir="wasm_outputs"
mkdir -p "$output_dir"

# Iterate over each tool and compile to WebAssembly
for tool in "${tools[@]}"; do
    tool_dir="src/$tool"
    cpp_file="$tool_dir/${tool}.cpp"
    wasm_output="$output_dir/${tool}.html"
    
    echo "Compiling $tool to WebAssembly..."
    
    emcc "$cpp_file" -o "$wasm_output" \
        -s WASM=1 \
        -s "EXPORTED_FUNCTIONS=['_main']" \
        --no-entry \
        -O3

    echo "$tool compiled successfully: $wasm_output"
done

echo "All tools have been compiled to WebAssembly and are available in the '$output_dir' directory."
