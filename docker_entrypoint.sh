#!/usr/bin/env bash
# Entrypoint for VCFX Docker image.
# It adds VCFX tool directories to the PATH and then executes the given command.

# Source the helper script if available
if [ -f /usr/local/bin/add_vcfx_tools_to_path.sh ]; then
    source /usr/local/bin/add_vcfx_tools_to_path.sh
fi

exec "$@"
