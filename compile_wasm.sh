#!/usr/bin/env bash
set -e

mkdir -p build_wasm
cd build_wasm

# Turn on BUILD_WASM
cmake -DBUILD_WASM=ON ..

cmake --build .

echo "All VCFX tools built for WebAssembly in build_wasm/."
echo "Use 'ls -R build_wasm' to see output. If you want .html or .js from Emscripten, you can adjust linking flags or suffixes."
