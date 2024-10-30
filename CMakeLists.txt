cmake_minimum_required(VERSION 3.10)
project(VCFX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Help find zlib on macOS with Homebrew
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/zlib")
    # For Apple Silicon (M1/M2)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/zlib")
endif()

find_package(ZLIB REQUIRED)

# ARM-specific optimizations
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -O3")
endif()

# Add clean target
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_clean.cmake
)

include_directories(include)

add_subdirectory(src)
add_subdirectory(tests)
