#!/usr/bin/env bash
#
# Adds all built VCFX tool directories in build/src to the PATH, so you can
# invoke them by name.
#
# Usage:
#   source ./add_vcfx_tools_to_path.sh

# Where is the root of this script? (i.e., your VCFX repository root)
# Adjust if needed; for example if you keep this script in the top-level dir:
REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Our compiled tools should be under build/src
BUILD_SRC_DIR="${REPO_ROOT}/build/src"

# Check that this path exists:
if [ ! -d "${BUILD_SRC_DIR}" ]; then
    echo "Error: build/src directory not found at: ${BUILD_SRC_DIR}"
    echo "Make sure you have run 'cmake .. && make' inside ./build"
    return 1
fi

# We'll gather a list of directories under build/src/VCFX_*
# that actually contain an executable matching the pattern "VCFX_*"
# Then add those directories to PATH.

TOOL_DIRS=""
while IFS= read -r -d '' toolExec; do
    # 'toolExec' is something like: build/src/VCFX_af_subsetter/VCFX_af_subsetter
    toolDir=$(dirname "$toolExec")
    # Only add it once if not present
    if [[ ":$TOOL_DIRS:" != *":$toolDir:"* ]]; then
        TOOL_DIRS="${TOOL_DIRS}:${toolDir}"
    fi
done < <(find "${BUILD_SRC_DIR}" -type f -perm /111 -name 'VCFX_*' -print0 2>/dev/null)

# If empty (no tools found), bail out
if [ -z "$TOOL_DIRS" ]; then
    echo "Warning: No VCFX tools found in ${BUILD_SRC_DIR}. Did you run 'make'?"
else
    # Remove leading colon
    TOOL_DIRS="${TOOL_DIRS#:}"

    # Add them to PATH
    export PATH="${TOOL_DIRS}:${PATH}"

    echo "Added these VCFX tool directories to PATH:"
    echo "${TOOL_DIRS//:/\\n}"
fi

# End of script
