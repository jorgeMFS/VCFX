# tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.14)

# We can build a small test or multiple tests
# For example, a simple program test_basic or we can do direct add_test calls

add_executable(test_basic test_basic.cpp)
target_link_libraries(test_basic PRIVATE vcfx_core)

# We define some ctest-based test:
add_test(NAME BasicTest COMMAND test_basic)

# We can also do direct testing of tools:
# e.g. check that 'VCFX_variant_counter' runs with a small VCF and returns 0
add_test(NAME variant_counter_runs
    COMMAND VCFX_variant_counter
    < ${CMAKE_CURRENT_SOURCE_DIR}/test_data/input_small.vcf
)
