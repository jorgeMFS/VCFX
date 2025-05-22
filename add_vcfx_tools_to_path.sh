#!/usr/bin/env bash
#
# Adds all built VCFX tool directories in build/src to the PATH, so you can
# invoke them by name.
#
# Usage:
#   source ./add_vcfx_tools_to_path.sh

# Determine potential base directories that may contain VCFX tools.
# When running from the build tree this will be build/src, but inside the
# Docker image the tools reside in /usr/local/bin/VCFX_*/.
REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

BASE_DIRS=()
BUILD_SRC_DIR="${REPO_ROOT}/build/src"
if [ -d "${BUILD_SRC_DIR}" ]; then
    BASE_DIRS+=("${BUILD_SRC_DIR}")
fi

# Also check the standard installation prefix used in the Docker image
if compgen -G "/usr/local/bin/VCFX_*" > /dev/null; then
    BASE_DIRS+=("/usr/local/bin")
fi

if [ ${#BASE_DIRS[@]} -eq 0 ]; then
    echo "Warning: No VCFX tool directories found."
    return 1
fi

# Gather directories containing executables named VCFX_*
TOOL_DIRS=""
for base in "${BASE_DIRS[@]}"; do
    while IFS= read -r -d '' toolExec; do
        toolDir=$(dirname "$toolExec")
        if [[ ":$TOOL_DIRS:" != *":$toolDir:"* ]]; then
            TOOL_DIRS="${TOOL_DIRS}:${toolDir}"
        fi
    done < <(find "$base" -type f -perm /111 -name 'VCFX_*' -print0 2>/dev/null)
done

# If empty (no tools found), bail out
if [ -z "$TOOL_DIRS" ]; then
    echo "Warning: No VCFX tools found."
else
    # Remove leading colon
    TOOL_DIRS="${TOOL_DIRS#:}"

    # Add them to PATH
    export PATH="${TOOL_DIRS}:${PATH}"

    echo "Added these VCFX tool directories to PATH:"
    echo "${TOOL_DIRS//:/\\n}"
fi

# End of script
