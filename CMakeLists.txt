cmake_minimum_required(VERSION 3.14)
project(VCFX CXX)

# Option to build WebAssembly with Emscripten
option(BUILD_WASM "Build with emscripten toolchain" OFF)

if(BUILD_WASM)
    # Adjust this path for your Emscripten installation:
    set(CMAKE_TOOLCHAIN_FILE "/path/to/emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake" CACHE FILEPATH "Emscripten toolchain" FORCE)
    message(STATUS "Building for WebAssembly using Emscripten toolchain.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# If needed, find zlib or other libs
# find_package(ZLIB REQUIRED)

# If building natively on ARM
if((NOT BUILD_WASM) AND (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -O3")
endif()

enable_testing()

# We keep a top-level 'src' folder that has all subfolders for each tool
add_subdirectory(src)

# A separate directory for unit tests
add_subdirectory(tests)
