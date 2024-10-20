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
