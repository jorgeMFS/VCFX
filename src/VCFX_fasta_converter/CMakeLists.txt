add_executable(VCFX_fasta_converter VCFX_fasta_converter.cpp VCFX_fasta_converter.h)
target_link_libraries(VCFX_fasta_converter PRIVATE vcfx_core)
