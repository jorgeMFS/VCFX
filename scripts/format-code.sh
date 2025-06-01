#!/bin/bash
# Format all C/C++ source files using clang-format

set -e

echo "Running clang-format on C/C++ source files..."

# Find all C/C++ source files (excluding build directory)
find . -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" | \
    grep -v build | \
    xargs clang-format -i --style=file

echo "Code formatting complete."

# Check if any files were modified
if [[ -n $(git status --porcelain) ]]; then
    echo ""
    echo "The following files were formatted:"
    git status --porcelain
    echo ""
    echo "Please review the changes and commit them if they look correct:"
    echo "  git add -A"
    echo "  git commit -m \"Apply clang-format to C/C++ source files\""
else
    echo "No files needed formatting."
fi 