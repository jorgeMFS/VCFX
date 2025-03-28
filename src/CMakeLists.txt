cmake_minimum_required(VERSION 3.14)

# Build the core library from your shared code
add_library(vcfx_core STATIC vcfx_core.cpp)
target_include_directories(vcfx_core PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../include)

# Add all tool subdirectories
add_subdirectory(VCFX_header_parser)
add_subdirectory(VCFX_record_filter)
add_subdirectory(VCFX_field_extractor)
add_subdirectory(VCFX_format_converter)
add_subdirectory(VCFX_variant_counter)
add_subdirectory(VCFX_sample_extractor)
add_subdirectory(VCFX_sorter)
add_subdirectory(VCFX_validator)
add_subdirectory(VCFX_subsampler)
add_subdirectory(VCFX_genotype_query)
add_subdirectory(VCFX_allele_freq_calc)
add_subdirectory(VCFX_indexer)
add_subdirectory(VCFX_compressor)
add_subdirectory(VCFX_position_subsetter)
add_subdirectory(VCFX_haplotype_extractor)
add_subdirectory(VCFX_info_parser)
add_subdirectory(VCFX_variant_classifier)
add_subdirectory(VCFX_duplicate_remover)
add_subdirectory(VCFX_info_summarizer)
add_subdirectory(VCFX_distance_calculator)
add_subdirectory(VCFX_multiallelic_splitter)
add_subdirectory(VCFX_missing_data_handler)
add_subdirectory(VCFX_concordance_checker)
add_subdirectory(VCFX_allele_balance_calc)
add_subdirectory(VCFX_allele_counter)
add_subdirectory(VCFX_phase_checker)
add_subdirectory(VCFX_annotation_extractor)
add_subdirectory(VCFX_phred_filter)
add_subdirectory(VCFX_merger)
add_subdirectory(VCFX_metadata_summarizer)
add_subdirectory(VCFX_hwe_tester)
add_subdirectory(VCFX_fasta_converter)
add_subdirectory(VCFX_nonref_filter)
add_subdirectory(VCFX_dosage_calculator)
add_subdirectory(VCFX_population_filter)
add_subdirectory(VCFX_file_splitter)
add_subdirectory(VCFX_gl_filter)
add_subdirectory(VCFX_ref_comparator)
add_subdirectory(VCFX_ancestry_inferrer)
add_subdirectory(VCFX_impact_filter)
add_subdirectory(VCFX_info_aggregator)
add_subdirectory(VCFX_probability_filter)
add_subdirectory(VCFX_diff_tool)
add_subdirectory(VCFX_cross_sample_concordance)
add_subdirectory(VCFX_phase_quality_filter)
add_subdirectory(VCFX_indel_normalizer)
add_subdirectory(VCFX_custom_annotator)
add_subdirectory(VCFX_region_subsampler)
add_subdirectory(VCFX_allele_balance_filter)
add_subdirectory(VCFX_missing_detector)
add_subdirectory(VCFX_haplotype_phaser)
add_subdirectory(VCFX_af_subsetter)
add_subdirectory(VCFX_sv_handler)
add_subdirectory(VCFX_reformatter)
add_subdirectory(VCFX_quality_adjuster)
add_subdirectory(VCFX_inbreeding_calculator)
add_subdirectory(VCFX_outlier_detector)
add_subdirectory(VCFX_alignment_checker)
add_subdirectory(VCFX_ancestry_assigner)
add_subdirectory(VCFX_ld_calculator)

# Install the core library
install(TARGETS vcfx_core
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Define a list of all tool executables for installation
set(VCFX_TOOLS
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

# Install all tool executables
install(TARGETS ${VCFX_TOOLS}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

message(STATUS "VCFX tools to be installed: ${VCFX_TOOLS}")