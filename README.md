# VCFX: Comprehensive VCF Manipulation Toolkit

<p align="center">
  <img src="assets/images/VCFX.png" alt="VCFX Logo" width="180"/>
</p>

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://ieeta-pt.github.io/VCFX/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/vcfx.svg)](https://anaconda.org/bioconda/vcfx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker](https://img.shields.io/badge/Docker-GHCR-blue)](https://ieeta-pt.github.io/VCFX/docker/)
[![PyPI version](https://img.shields.io/pypi/v/vcfx.svg)](https://pypi.org/project/vcfx/)
[![PyPI downloads](https://img.shields.io/pypi/dm/vcfx.svg)](https://pypi.org/project/vcfx/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vcfx.svg)](https://pypi.org/project/vcfx/)

VCFX is a set of small C/C++ command line tools for manipulating and analysing Variant Call Format (VCF) files. Each tool does one job well and they can be chained together using standard streams.

## Features

- **60+ Specialized Tools** for filtering, transforming and analysing variants
- **Pipeline Ready**: tools read from stdin and write to stdout
- **Fast**: designed for large genomic datasets
- **Cross Platform**: Linux and macOS support
- **WebAssembly Builds** for browser or Node.js usage
- **Easy Installation** via PyPI, Bioconda or Docker
- **Python Bindings** for programmatic access

## Installation

### PyPI (Python Package)
```bash
pip install vcfx
```
After installing the Python bindings you can run any tool directly:
```python
import vcfx
vcfx.run_tool("alignment_checker", "--help")
```

### Bioconda
```bash
conda install -c bioconda vcfx
```

### Docker
```bash
docker pull ghcr.io/jorgemfs/vcfx:latest
docker run --rm ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name --help
```

### Build from Source
```bash
git clone https://github.com/ieeta-pt/VCFX.git
cd VCFX
mkdir build && cd build
cmake .. -DPYTHON_BINDINGS=ON
make
```
Optionally run `make install` to place the tools in `~/.local/bin`.

## Quick Example
```bash
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_allele_freq_calc > snp_frequencies.tsv
```

## Documentation

Full documentation is available at [ieeta-pt.github.io/VCFX](https://ieeta-pt.github.io/VCFX/). The `docs` folder contains the sources and can be served locally with `mkdocs serve`.

## Repository Layout

- `src/` – C++ source code for all tools
- `python/` – optional Python bindings
- `docs/` – documentation sources
- `tests/` – shell and Python tests
- `examples/` – usage examples and workflows

## Development

Run tests with `pytest` and `ctest` from the build directory. Code style is enforced with `clang-format` and pre-commit hooks:
```bash
pre-commit install
```

## Benchmarks

VCFX tools are optimized for high-throughput genomic analysis. Below are benchmark results on a **4.0 GB VCF file** (chromosome 21 from 1000 Genomes Phase 3, ~427K variants, 2,504 samples).

### Performance Summary by Category

| Category | Tools | Avg Time (s) | Min (s) | Max (s) |
|----------|------:|-------------:|--------:|--------:|
| Basic I/O & Validation | 4 | 2.2 | 0.10 | 8.5 |
| Filtering | 8 | 7.6 | 0.10 | 14.4 |
| Transformation | 7 | 6.6 | 0.09 | 20.0 |
| Extraction & Subsetting | 6 | 10.6 | 0.11 | 14.1 |
| Quality Control | 5 | 18.7 | 0.13 | 36.6 |
| Annotation & INFO | 5 | 28.1 | 9.7 | 35.1 |
| Analysis & Calculation | 8 | 40.4 | 8.2 | 147.3 |
| Haplotype | 2 | 39.4 | 25.1 | 53.7 |
| Ancestry | 2 | 20.3 | 0.11 | 40.6 |
| File Management | 5 | 42.5 | 0.09 | 186.3 |
| Format Conversion | 2 | 28.1 | 17.0 | 39.2 |
| Comparison & Diff | 3 | 76.1 | 15.4 | 187.9 |
| Advanced Analysis | 2 | 629.7 | 0.14 | 1259.2 |

### Comparison with Standard Tools

| Tool | Operation | Time (s) |
|------|-----------|----------|
| **VCFX** variant_counter | Count variants | 8.5 |
| **VCFX** allele_freq_calc | Allele frequencies | 29.1 |
| **VCFX** missing_detector | Missing data | 9.3 |
| bcftools view | Parse VCF | 37.9 |
| bcftools query | Extract fields | 2.2 |
| bcftools filter | Filter variants | 19.4 |
| vcftools --freq | Allele frequencies | 86.2 |
| vcftools --missing | Missing data | 82.9 |

VCFX achieves competitive or better performance than bcftools/vcftools for many common operations, with the advantage of modular single-purpose tools.

### Running Benchmarks

```bash
cd benchmarks
make datasets
make tools-check
make all
```

Results are saved to `benchmarks/results/` and visualization notebooks are under `benchmarks/notebooks/`.

## Citation

If you use VCFX in your research please cite:
```
@inproceedings{silva2025vcfx,
  title={VCFX: A Minimalist, Modular Toolkit for Streamlined Variant Analysis},
  author={Silva, Jorge Miguel and Oliveira, José Luis},
  booktitle={12th International Work-Conference on Bioinformatics and Biomedical Engineering (IWBBIO 2025)},
  year={2025},
  organization={Springer}
}
```

## License

VCFX is distributed under the [MIT License](LICENSE).
