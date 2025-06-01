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

## Repository Structure

- `src/` - C++ source code for all command-line tools.
- `docs/` - MkDocs documentation sources.
- `python/` - Optional Python bindings and helper modules.
- `tests/` - Shell and Python tests with sample data.
- `examples/` - Small usage examples and workflow demos.
- `scripts/` - Utility scripts for development and CI.
- `include/` - Public C++ headers.
- `assets/` - Project logos and other media assets.

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
cmake .. -DPYTHON_BINDINGS=ON
make
```

Enabling `PYTHON_BINDINGS` builds the Python module.

### Versioning

The version numbers for both the C++ tools and the Python package are defined
once in `CMakeLists.txt` (`VCFX_VERSION_MAJOR`, `VCFX_VERSION_MINOR`,
`VCFX_VERSION_PATCH`). The Python build reads these values at install time to
ensure a single source of truth.
For a detailed list of changes in each release, see [CHANGELOG.md](CHANGELOG.md).

#### Installing Python Bindings from a Local Checkout

```bash
python3 -m venv venv && source venv/bin/activate
pip install --no-build-isolation ./python
```

When offline, ensure the `setuptools` and `cmake` packages are already installed
because build isolation is disabled.

#### Building a Wheel

To build a distributable wheel from the local checkout:

```bash
python3 -m venv venv && source venv/bin/activate
python3 -m pip wheel ./python -w dist --no-build-isolation
```

When offline, install `setuptools` and `wheel` beforehand:

```bash
python3 -m pip install setuptools wheel
```

#### Installing from PyPI or a Wheel

```bash
pip install vcfx               # from PyPI
pip install dist/vcfx-*.whl    # from a local wheel
```

### Python API

VCFX also offers optional Python bindings. See
[the Python API docs](docs/python_api.md) for details. After installing, you can
list available command line tools via `vcfx.available_tools()` and run any tool
with `vcfx.run_tool()`:

```python
import vcfx
vcfx.run_tool("alignment_checker", "--help")
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
### Common Flags

All tools recognize `-h`/`--help` and `-v`/`--version` via `vcfx::handle_common_flags`.
A minimal `main()` typically looks like:
```cpp
static void show_help() { printHelp(); }

int main(int argc, char* argv[]) {
    if (vcfx::handle_common_flags(argc, argv, "VCFX_example", show_help))
        return 0;
    // ...
}
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

Python tests are executed with `pytest`. Install the development requirements first:

```bash
pip install -r requirements-dev.txt
```

After building the bindings using `cmake -DPYTHON_BINDINGS=ON` (or after installing the wheel), run:

```bash
pytest tests/python
```

To automatically run `ruff`, `flake8`, `mypy`, and `pytest` before each commit,
install the pre-commit hooks once:

```bash
pre-commit install
```


## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/CONTRIBUTING.md) for details.

## Citation

If you use VCFX in your research, please cite (see also [CITATION.cff](CITATION.cff)):

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

## Development

### Code Formatting

This project uses `clang-format` for C/C++ code formatting. The clang-format hook is skipped in GitHub Actions CI to avoid version conflicts, but developers should run it locally before committing.

#### Installing clang-format

**macOS:**
```bash
brew install clang-format
```

**Ubuntu/Debian:**
```bash
sudo apt-get install clang-format
```

#### Formatting Code

To format all C/C++ source files:
```bash
# Using the provided script
./scripts/format-code.sh

# Or manually
find . -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" | \
  grep -v build | \
  xargs clang-format -i --style=file
```

#### Pre-commit Hooks

The project uses pre-commit hooks for code quality. Install them with:
```bash
pre-commit install
```

**Note:** If clang-format is not available on your system, you can skip it:
```bash
SKIP=clang-format pre-commit run --all-files
```
