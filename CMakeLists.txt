cmake_minimum_required(VERSION 3.14)

# Set version information
set(VCFX_VERSION_MAJOR 1)
set(VCFX_VERSION_MINOR 0)
set(VCFX_VERSION_PATCH 2)
set(VCFX_VERSION "${VCFX_VERSION_MAJOR}.${VCFX_VERSION_MINOR}.${VCFX_VERSION_PATCH}")

project(VCFX 
  VERSION ${VCFX_VERSION}
  DESCRIPTION "A Comprehensive VCF Manipulation Toolkit"
  LANGUAGES CXX
)

# Optionally allow building for WebAssembly via Emscripten
option(BUILD_WASM "Build with emscripten toolchain" OFF)

if(BUILD_WASM)
    set(CMAKE_TOOLCHAIN_FILE "/path/to/emscripten.cmake" CACHE FILEPATH "Emscripten toolchain" FORCE)
    message(STATUS "Building for WebAssembly (Emscripten).")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Locate zlib if we need it (some tools do)
# On macOS with Homebrew, you might need to specify paths. E.g.:
# list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/zlib" "/opt/homebrew/opt/zlib")
find_package(ZLIB REQUIRED)

# If building natively on ARM, e.g. Apple M1
if((NOT BUILD_WASM) AND (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -O3")
endif()

enable_testing()

# Add top-level 'src' subdirectory, which in turn references each tool subdirectory
add_subdirectory(src)

# Add a tests subdir if you have tests
# Comment out this line since we don't have a CMakeLists.txt file in the tests directory
# add_subdirectory(tests)

# Installation configuration
include(GNUInstallDirs)

# Install header files
install(FILES 
  include/vcfx_core.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/vcfx
)

# Install documentation files if they exist
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
  install(FILES 
    README.md
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
  )
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
  install(FILES 
    LICENSE
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
  )
endif()

message(STATUS "VCFX version: ${VCFX_VERSION}")
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
