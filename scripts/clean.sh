#!/bin/bash

# Remove build directory
rm -rf build/*

# Remove any temporary files
find . -name "*.o" -type f -delete
find . -name "*.a" -type f -delete
find . -name "*.so" -type f -delete
find . -name "*.dylib" -type f -delete
find . -name "*.dll" -type f -delete
find . -name "CMakeCache.txt" -type f -delete
find . -name "CMakeFiles" -type d -exec rm -rf {} +
find . -name "Makefile" -type f -delete
find . -name "cmake_install.cmake" -type f -delete

# Remove WebAssembly outputs
rm -rf wasm_outputs/

echo "Clean completed successfully!"