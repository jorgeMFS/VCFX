# VCFX: Comprehensive VCF Manipulation Toolkit

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://jorgeMFS.github.io/VCFX/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/vcfx.svg)](https://anaconda.org/bioconda/vcfx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

VCFX is a collection of specialized C/C++ command-line tools designed for efficient manipulation, analysis, and transformation of VCF (Variant Call Format) files used in genomic research and bioinformatics. Each tool is an independent executable that follows the Unix philosophy: do one thing well and work seamlessly with other tools through standard input/output streams.

## Documentation

**[View the full documentation](https://jorgeMFS.github.io/VCFX/)**

For local documentation, see the `docs` directory or build the documentation site:

```bash
pip install mkdocs-material pymdown-extensions
mkdocs serve
```

## Key Features

- **60 Specialized Tools**: Each optimized for a specific VCF-related task
- **Pipeline-Ready**: All tools work with standard input/output for easy integration
- **Performance-Focused**: Designed for handling large genomic datasets efficiently
- **Cross-Platform**: Works on Linux, macOS, and Windows
- **WebAssembly Support**: Optional WASM builds for browser/Node.js environments
- **Bioconda Available**: Easy installation through Bioconda package manager

## Tool Categories

- **Data Analysis**: Extract statistical information from variant data
- **Data Filtering**: Select variants based on specific criteria
- **Data Transformation**: Convert or reformat VCF data
- **Quality Control**: Validate and check data quality
- **File Management**: Handle VCF files efficiently
- **Annotation and Reporting**: Add or extract annotations from VCF files
- **Data Processing**: Process variants and samples

## Quick Start

### Installation

#### Using Bioconda

```bash
conda install -c bioconda vcfx
```

#### Building from Source

```bash
git clone https://github.com/jorgeMFS/VCFX.git
cd VCFX
mkdir -p build && cd build
cmake ..
make
```

### Basic Usage Example

```bash
# Calculate allele frequencies for SNPs only
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_allele_freq_calc > snp_frequencies.tsv
```

## Building for WebAssembly

If you have [Emscripten](https://emscripten.org/) installed:

```bash
mkdir build_wasm && cd build_wasm
cmake -DBUILD_WASM=ON ..
cmake --build .
```

## Running Tests

```bash
cd build
ctest --verbose
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.