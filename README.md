<p align="center">
  <img src="assets/images/VCFX.png" alt="VCFX Logo" width="180"/>
</p>

# Comprehensive VCF Manipulation Toolkit


[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://ieeta-pt.github.io/VCFX/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/vcfx.svg)](https://anaconda.org/bioconda/vcfx)
[![Platforms](https://anaconda.org/bioconda/vcfx/badges/platforms.svg)](https://anaconda.org/bioconda/vcfx)
[![Downloads](https://anaconda.org/bioconda/vcfx/badges/downloads.svg)](https://anaconda.org/bioconda/vcfx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker](https://img.shields.io/badge/Docker-GHCR-blue)](https://ieeta-pt.github.io/VCFX/docker/)

VCFX is a collection of specialized C/C++ command-line tools designed for efficient manipulation, analysis, and transformation of VCF (Variant Call Format) files used in genomic research and bioinformatics. Each tool is an independent executable that follows the Unix philosophy: do one thing well and work seamlessly with other tools through standard input/output streams.

## Documentation

**[View the full documentation](https://ieeta-pt.github.io/VCFX/)**

For local documentation, see the `docs` directory or build the documentation site:

```bash
pip install mkdocs-material pymdown-extensions
mkdocs serve
```

## Key Features

- **60 Specialized Tools**: Each optimized for a specific VCF-related task
- **Pipeline-Ready**: All tools work with standard input/output for easy integration
- **Performance-Focused**: Designed for handling large genomic datasets efficiently
- **Cross-Platform**: Works on Linux and macOS
- **WebAssembly Support**: Optional WASM builds for browser/Node.js environments
- **Bioconda Available**: Easy installation through Bioconda package manager
- **Docker Support**: Ready-to-use Docker image available via GitHub Container Registry

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
git clone https://github.com/ieeta-pt/VCFX.git
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

### Listing Available Tools

```bash
vcfx list
```

### Show Tool Documentation

```bash
vcfx help allele_counter
```

## Building for WebAssembly

If you have [Emscripten](https://emscripten.org/) installed:

```bash
mkdir build_wasm && cd build_wasm
emcmake cmake -DBUILD_WASM=ON ..
cmake --build .
```

## Running Tests

From your build directory, run:

```bash
ctest --output-on-failure
```

You can also execute all shell scripts directly with:

```bash
bash ../tests/test_all.sh
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/CONTRIBUTING.md) for details.

## Citation

If you use VCFX in your research, please cite:

```
@inproceedings{silva2025vcfx,
  title={VCFX: A Minimalist, Modular Toolkit for Streamlined Variant Analysis},
  author={Silva, Jorge Miguel and Oliveira, Jos{\'e} Luis},
  booktitle={12th International Work-Conference on Bioinformatics and Biomedical Engineering (IWBBIO 2025)},
  year={2025},
  organization={Springer}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
