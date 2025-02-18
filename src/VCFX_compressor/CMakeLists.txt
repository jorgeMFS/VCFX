cmake_minimum_required(VERSION 3.14)
project(VCFX_compressor LANGUAGES CXX)

add_executable(VCFX_compressor
    VCFX_compressor.cpp
    VCFX_compressor.h
)

# Link against the core library + ZLIB
target_link_libraries(VCFX_compressor PRIVATE vcfx_core ZLIB::ZLIB)
