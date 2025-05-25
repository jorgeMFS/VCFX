#!/usr/bin/env bash
set -e

mkdir -p build_wasm
cd build_wasm

# Turn on BUILD_WASM using emcmake if available
if command -v emcmake >/dev/null 2>&1; then
    emcmake cmake -DBUILD_WASM=ON ..
else
    cmake -DBUILD_WASM=ON ..
fi

cmake --build .

echo "All VCFX tools and the vcfx wrapper built for WebAssembly in build_wasm/."
echo "Use 'ls -R build_wasm' to see output. If you want .html or .js from Emscripten, you can adjust linking flags or suffixes."
