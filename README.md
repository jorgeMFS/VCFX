# VCFX: Comprehensive VCF Manipulation Toolkit

<p align="center">
  <img src="assets/images/VCFX.png" alt="VCFX Logo" width="400"/>
</p>

<p align="center">
  <a href="https://github.com/jorgeMFS/VCFX/actions"><img src="https://github.com/jorgeMFS/VCFX/actions/workflows/docs.yml/badge.svg" alt="Documentation Status"></a>
  <a href="https://anaconda.org/bioconda/vcfx"><img src="https://img.shields.io/conda/vn/bioconda/vcfx.svg" alt="Conda Version"></a>
  <a href="https://github.com/jorgeMFS/VCFX/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License"></a>
</p>

VCFX is a collection of specialized C/C++ command-line tools designed for efficient manipulation, analysis, and transformation of VCF (Variant Call Format) files used in genomic research and bioinformatics. Each tool is an independent executable that follows the Unix philosophy: do one thing well and work seamlessly with other tools through standard input/output streams.

## Documentation

**[View the full documentation](https://jorgeMFS.github.io/VCFX/)**

For local documentation, see the `docs` directory or build the documentation site:
```bash
pip install mkdocs-material pymdown-extensions
mkdocs serve
```

## Installation

### Using Bioconda (Recommended)

```bash
# Set up Bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install VCFX
conda install vcfx
```

### Building from Source

```bash
git clone https://github.com/jorgeMFS/VCFX.git
cd VCFX
mkdir -p build && cd build
cmake ..
make
```

See the [Installation Guide](https://jorgeMFS.github.io/VCFX/installation) for more options.

## Key Features

- **60 Specialized Tools**: Each optimized for a specific VCF-related task
- **Pipeline-Ready**: All tools work with standard input/output for easy integration
- **Performance-Focused**: Designed for handling large genomic datasets efficiently
- **Cross-Platform**: Works on Linux, macOS, and Windows
- **WebAssembly Support**: Optional WASM builds for browser/Node.js environments
- **Bioconda Package**: Easy installation via conda

## Tool Categories

- **Data Analysis**: Extract statistical information from variant data
- **Data Filtering**: Select variants based on specific criteria
- **Data Transformation**: Convert or reformat VCF data
- **Quality Control**: Validate and check data quality
- **File Management**: Handle VCF files efficiently
- **Annotation and Reporting**: Add or extract annotations from VCF files
- **Data Processing**: Process variants and samples

## Quick Start Example

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

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/CONTRIBUTING.md) for details.

## Running Tests

VCFX includes a comprehensive test suite with dedicated test scripts for each tool. You can run either all tests or individual tool tests.

### Running All Tests

```bash
# Navigate to the tests directory
cd tests

# Make the test script executable if needed
chmod +x test_all.sh

# Run all tests
./test_all.sh
```

### Running Individual Tool Tests

To test a specific tool:

```bash
# Navigate to the tests directory
cd tests

# Make the test script executable
chmod +x test_tool_name.sh

# Run the specific test
./test_tool_name.sh
```

Replace `test_tool_name.sh` with the appropriate test script (e.g., `test_allele_freq_calc.sh`).

### Test Structure

The test suite uses:
- Sample data in the `tests/data` directory
- Expected outputs in the `tests/expected` directory
- Test output is saved to the `tests/out` directory

Each test script validates that the tool produces the expected output for various input scenarios.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.