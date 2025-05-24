# Python API

VCFX provides optional Python bindings exposing a subset of helper
functions from the C++ `vcfx_core` library. The bindings are built as a
native Python extension and can be enabled through CMake.

## Installation

Build the project with the `PYTHON_BINDINGS` option enabled:

```bash
mkdir build && cd build
cmake -DPYTHON_BINDINGS=ON ..
make -j
```

The compiled module will be placed in the `build/python` directory.
You can also install the package via `pip` which will invoke CMake
automatically:

```bash
pip install ./python
```

## Available Functions

The module exposes the following helpers:

- `trim(text)` – remove leading and trailing whitespace.
- `split(text, delimiter)` – split `text` by the given delimiter and
  return a list of strings.
- `read_file_maybe_compressed(path)` – read a plain or gzip/BGZF
  compressed file and return its contents as a string.
- `get_version()` – return the VCFX version string.

## Example Usage

```python
import vcfx

print(vcfx.trim("  abc  "))
# 'abc'

print(vcfx.split("A,B,C", ","))
# ['A', 'B', 'C']

version = vcfx.get_version()
print("VCFX version:", version)
```

## Tool Wrappers

Besides the helper functions, the package provides lightweight wrappers for
all command line tools shipped with VCFX.  The wrappers simply invoke the
corresponding ``VCFX_*`` executable via ``subprocess``.

Use ``vcfx.available_tools()`` to see which tools are accessible on your
``PATH`` and call them either via ``vcfx.run_tool(name, *args)`` or by using
the tool name as a function:

```python
import vcfx

print(vcfx.available_tools())

# run through the generic helper
vcfx.run_tool("alignment_checker", "--help")

# or directly by name (if available)
vcfx.alignment_checker("--help")
```
