#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
BUILD_DIR="${ROOT_DIR}/build/python_tool_wrappers"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -DPYTHON_BINDINGS=ON ../..
make -j

# ensure tools on PATH for the wrappers
for exe in "${BUILD_DIR}"/src/VCFX_*/*; do
    if [ -x "$exe" ]; then
        PATH="$(dirname "$exe"):${PATH}"
    fi
done
export PATH

cd "$SCRIPT_DIR"

PYTHONPATH="${BUILD_DIR}/python" python3 - <<'PY'
import vcfx

rows = vcfx.alignment_checker("data/align_Y.vcf", "data/align_refY.fa")
assert rows and rows[0]["Discrepancy_Type"] == "ALT_MISMATCH"

counts = vcfx.allele_counter("data/allele_counter_A.vcf")
assert counts[0]["Sample"] == "S1"

n = vcfx.variant_counter("data/variant_counter_normal.vcf")
assert n == 5

freqs = vcfx.allele_freq_calc("data/allele_freq_calc/simple.vcf")
assert freqs[0]["Allele_Frequency"] == "0.5000"

annotated = vcfx.info_aggregator("data/aggregator/basic.vcf", ["DP"])
assert "#AGGREGATION_SUMMARY" in annotated

parsed = vcfx.info_parser("data/info_parser/basic.vcf", ["DP"])
assert parsed[0]["DP"] == "10"

summary = vcfx.info_summarizer("data/info_summarizer/basic.vcf", ["DP"])
assert summary[0]["Mean"] == "20.0000"

fasta = vcfx.fasta_converter("data/fasta_converter/basic.vcf")
assert fasta.startswith(">")
print("Python tool wrappers OK")
PY
