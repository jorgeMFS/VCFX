cmake_minimum_required(VERSION 3.14)
project(VCFX_variant_counter LANGUAGES CXX)

add_executable(VCFX_variant_counter
    VCFX_variant_counter.cpp
    VCFX_variant_counter.h
)

target_link_libraries(VCFX_variant_counter PRIVATE vcfx_core)
