#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define an array of tool names
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

# Create main directories
mkdir -p src include tests docs examples scripts

# Create src subdirectories for each tool and corresponding files
for tool in "${tools[@]}"; do
    tool_dir="src/$tool"
    mkdir -p "$tool_dir"
    
    # Create CMakeLists.txt for each tool
    cat << EOF > "$tool_dir/CMakeLists.txt"
add_executable($tool ${tool}.cpp ${tool}.h)
target_link_libraries($tool PRIVATE vcfx_core)
EOF

    # Create an empty ${tool}.cpp file with a placeholder main function
    cat << EOF > "$tool_dir/${tool}.cpp"
#include "${tool}.h"

int main(int argc, char* argv[]) {
    // TODO: Implement $tool functionality
    return 0;
}
EOF

    # Create an empty ${tool}.h file with include guards
    guard=$(echo "$tool" | tr '[:lower:]' '[:upper:]')
    cat << EOF > "$tool_dir/${tool}.h"
#ifndef ${guard}_H
#define ${guard}_H

// Declarations for $tool

#endif // ${guard}_H
EOF
done

# Create main CMakeLists.txt
cat << EOF > CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(VCFX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_subdirectory(src)
add_subdirectory(tests)
EOF

# Create src/CMakeLists.txt with all tool subdirectories
cat << EOF > src/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

# Core library
add_library(vcfx_core src/vcfx_core.cpp)

# Add all tool subdirectories
EOF

for tool in "${tools[@]}"; do
    echo "add_subdirectory($tool)" >> src/CMakeLists.txt
done

# Create include/vcfx_core.h
cat << EOF > include/vcfx_core.h
#ifndef VCFX_CORE_H
#define VCFX_CORE_H

#include <iostream>
#include <string>

// Core functionalities for VCFX tools

#endif // VCFX_CORE_H
EOF

# Create src/vcfx_core.cpp
cat << EOF > src/vcfx_core.cpp
#include "vcfx_core.h"

// Implementation of core functionalities
EOF

# Create a main README.md
cat << EOF > README.md
# VCFX

VCFX is a collection of C/C++ tools for processing and analyzing VCF (Variant Call Format) files, with WebAssembly compatibility.

## Tools

$(for i in "${!tools[@]}"; do echo "$(($i + 1)). ${tools[$i]}"; done)

## Building

To build the project:

\`\`\`
mkdir build && cd build
cmake ..
make
\`\`\`

## Running Tests

To run the tests:

\`\`\`
cd build
ctest
\`\`\`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
EOF

# Create a LICENSE file (MIT License)
cat << EOF > LICENSE
MIT License

Copyright (c) $(date +%Y)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF

# Create a .gitignore file
cat << EOF > .gitignore
# Build directories
build/
bin/

# IDE-specific files
.vscode/
.idea/

# Compiled object files
*.o
*.obj

# Compiled dynamic libraries
*.so
*.dylib
*.dll

# Compiled static libraries
*.a
*.lib

# Executables
*.exe
*.out

# CMake generated files
CMakeCache.txt
CMakeFiles/
cmake_install.cmake
Makefile

# OS-specific files
.DS_Store
Thumbs.db

# Other
tools.md
prompt.md
EOF

echo "VCFX project structure with 60 tools has been set up successfully!"
