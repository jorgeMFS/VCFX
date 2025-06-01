#!/bin/bash
# Test PyPI package build locally

set -e

echo "Testing PyPI package build..."

# Extract version
VERSION=$(grep -E "project\(VCFX VERSION" CMakeLists.txt | sed -E 's/.*VERSION ([0-9.]+).*/\1/')
echo "Version from CMakeLists.txt: $VERSION"

# Create temporary build directory
TEMP_DIR=$(mktemp -d)
echo "Using temp directory: $TEMP_DIR"

# Copy Python package files
cp -r python/* "$TEMP_DIR/"
echo "$VERSION" > "$TEMP_DIR/VERSION"

# Build in temp directory
cd "$TEMP_DIR"
export VCFX_VERSION=$VERSION

echo "Building source distribution..."
python3 -m build --sdist

echo "Checking built package..."
ls -la dist/
python3 -m twine check dist/*

# Extract and check version from built package
echo "Verifying version in package..."
tar -xzf dist/*.tar.gz
PACKAGE_DIR=$(ls -d vcfx-*)
if [ -f "$PACKAGE_DIR/VERSION" ]; then
    PACKAGE_VERSION=$(cat "$PACKAGE_DIR/VERSION")
    echo "Version in package: $PACKAGE_VERSION"
    if [ "$PACKAGE_VERSION" = "$VERSION" ]; then
        echo "✅ Version matches!"
    else
        echo "❌ Version mismatch!"
        exit 1
    fi
else
    echo "❌ VERSION file not found in package!"
    exit 1
fi

# Cleanup
cd ..
rm -rf "$TEMP_DIR"

echo "✅ PyPI build test completed successfully!" 