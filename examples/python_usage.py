"""Example usage of the vcfx Python module.

This script shows how to source ``add_vcfx_tools_to_path.sh`` to make the
compiled tools available, query ``vcfx.available_tools()``, run a tool, and
perform basic error handling.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import vcfx

# ----------------------------------------------------------------------------
# 1. Source add_vcfx_tools_to_path.sh
# ----------------------------------------------------------------------------

# Determine the repository root relative to this file
REPO_ROOT = Path(__file__).resolve().parents[1]
script = REPO_ROOT / "add_vcfx_tools_to_path.sh"

# Use a subshell to source the script and print the resulting PATH
# so we can update our environment accordingly.
# When running interactively you can simply execute:
#   source /path/to/add_vcfx_tools_to_path.sh
cmd = ["bash", "-c", f"source {script} >/dev/null && printf '%s' \"$PATH\""]
new_path = subprocess.check_output(cmd, text=True)
os.environ["PATH"] = new_path

# ----------------------------------------------------------------------------
# 2. Query available tools and run one
# ----------------------------------------------------------------------------

print("Available tools:", vcfx.available_tools())

# Run a simple tool via the generic helper. Here we count the variants in a
# small VCF located under the tests directory.
vcf_file = REPO_ROOT / "tests" / "data" / "variant_counter_normal.vcf"

try:
    result = vcfx.run_tool("variant_counter", str(vcf_file), capture_output=True)
    print("variant_counter output:", result.stdout.strip())
except FileNotFoundError as exc:
    print("Tool not found:", exc)
except subprocess.CalledProcessError as exc:
    print("Tool failed with status", exc.returncode)


