#!/bin/bash

# Create main directories
mkdir -p src include tests docs examples scripts

# Create src subdirectories for each tool
for i in {1..60}; do
    mkdir -p "src/tool_$i"
done

# Create main CMakeLists.txt
cat << EOF > CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(VCFX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_subdirectory(src)
add_subdirectory(tests)
EOF

# Create src/CMakeLists.txt
cat << EOF > src/CMakeLists.txt
add_subdirectory(tool_1)
add_subdirectory(tool_2)
# Add more subdirectories as needed
EOF

# Create a sample tool CMakeLists.txt (for tool_1)
cat << EOF > src/tool_1/CMakeLists.txt
add_executable(tool_1 tool_1.cpp)
target_link_libraries(tool_1 PRIVATE vcfx_core)
EOF

# Create a sample header file
cat << EOF > include/vcfx_core.h
#ifndef VCFX_CORE_H
#define VCFX_CORE_H

// Core functionality for VCFX tools

#endif // VCFX_CORE_H
EOF

# Create a sample source file
cat << EOF > src/vcfx_core.cpp
#include "vcfx_core.h"

// Implementation of core functionality
EOF

# Create a main README.md
cat << EOF > README.md
# VCFX

VCFX is a collection of C/C++ tools for processing and analyzing VCF (Variant Call Format) files, with WebAssembly compatibility.

## Tools

1. VCF Header Parser
2. VCF Record Filter
3. VCF Field Extractor
...

## Building

To build the project:

\`\`\`
mkdir build && cd build
cmake ..
make
\`\`\`

## Running Tests

To run the tests:

\`\`\`
cd build
ctest
\`\`\`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
EOF

# Create a LICENSE file (MIT License)
cat << EOF > LICENSE
MIT License

Copyright (c) $(date +%Y) Your Name

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF

# Create a .gitignore file
cat << EOF > .gitignore
# Build directories
build/
bin/

# IDE-specific files
.vscode/
.idea/

# Compiled object files
*.o
*.obj

# Compiled dynamic libraries
*.so
*.dylib
*.dll

# Compiled static libraries
*.a
*.lib

# Executables
*.exe
*.out

# CMake generated files
CMakeCache.txt
CMakeFiles/
cmake_install.cmake
Makefile

# OS-specific files
.DS_Store
Thumbs.db

# Other
tools.md
prompt.md

EOF

# Initialize git repository
git init
git add .
git commit -m "Initial project setup"

echo "VCFX project structure has been set up successfully!"