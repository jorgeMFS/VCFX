# VCFX Python Package

Python bindings for the VCFX toolkit - a comprehensive VCF manipulation toolkit.

## Installation

```bash
pip install vcfx
```

## Quick Start

```python
import vcfx

# Use helper functions
text = vcfx.trim("  hello  ")  # Returns "hello"
parts = vcfx.split("A,B,C", ",")  # Returns ["A", "B", "C"]

# Get version
version = vcfx.get_version()
print(f"VCFX version: {version}")

# Use tool wrappers (requires VCFX tools in PATH)
count = vcfx.variant_counter("input.vcf")
freqs = vcfx.allele_freq_calc("input.vcf")
```

## Features

- **Native C++ bindings** for high-performance operations
- **Tool wrappers** for all VCFX command-line tools
- **Convenience functions** for common VCF analysis tasks
- **Type hints** for better development experience
- **Cross-platform support** (Linux, macOS)

## Requirements

- Python 3.10+
- For tool wrappers: VCFX command-line tools must be installed and available in PATH

## Documentation

For comprehensive documentation, visit: https://vcfx.readthedocs.io

## License

MIT License - see LICENSE file for details. 