# CMake tests for VCFX shell scripts

set(TEST_SCRIPTS
    test_af_subsetter.sh
    test_alignment_checker.sh
    test_allele_balance_calc.sh
    test_allele_balance_filter.sh
    test_allele_counter.sh
    test_allele_freq_calc.sh
    test_ancestry_assigner.sh
    test_ancestry_inferrer.sh
    test_annotation_extractor.sh
    test_compressor.sh
    test_concordance_checker.sh
    test_cross_sample_concordance.sh
    test_custom_annotator.sh
    test_diff_tool.sh
    test_distance_calculator.sh
    test_dosage_calculator.sh
    test_duplicate_remover.sh
    test_fasta_converter.sh
    test_field_extractor.sh
    test_file_splitter.sh
    test_format_converter.sh
    test_genotype_query.sh
    test_gl_filter.sh
    test_haplotype_extractor.sh
    test_header_parser.sh
    test_hwe_tester.sh
    test_impact_filter.sh
    test_indel_normalizer.sh
    test_indexer.sh
    test_info_aggregator.sh
    test_info_summarizer.sh
    test_inbreeding_calculator.sh
    test_ld_calculator.sh
    test_metadata_summarizer.sh
    test_merger.sh
    test_missing_data_handler.sh
    test_missing_detector.sh
    test_multiallelic_splitter.sh
    test_nonref_filter.sh
    test_outlier_detector.sh
    test_phase_checker.sh
    test_phase_quality_filter.sh
    test_phred_filter.sh
    test_population_filter.sh
    test_position_subsetter.sh
    test_probability_filter.sh
    test_quality_adjuster.sh
    test_record_filter.sh
    test_ref_comparator.sh
    test_reformatter.sh
    test_region_subsampler.sh
    test_sample_extractor.sh
    test_sorter.sh
    test_sv_handler.sh
    test_subsampler.sh
    test_validator.sh
    test_variant_classifier.sh
    test_variant_counter.sh
    test_doc_lookup.sh
    test_python_bindings.sh
    test_python_tool_wrappers.sh
    test_python_tool_wrappers_extra.sh
    test_python_version_fallback.sh
)

foreach(script ${TEST_SCRIPTS})
    get_filename_component(name ${script} NAME_WE)
    add_test(NAME ${name}
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/${script})
    set_tests_properties(${name} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()

# Optional Docker-based tests
find_program(DOCKER_EXECUTABLE docker)
if(DOCKER_EXECUTABLE)
    add_test(NAME test_docker
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/test_docker.sh)
    set_tests_properties(test_docker PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
