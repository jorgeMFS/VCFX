#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
BUILD_DIR="${ROOT_DIR}/build/python_tool_wrappers_extra"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -DPYTHON_BINDINGS=ON ../..
make -j

for exe in "${BUILD_DIR}"/src/VCFX_*/*; do
    if [ -x "$exe" ]; then
        PATH="$(dirname "$exe"):$PATH"
    fi
done
export PATH

cd "$SCRIPT_DIR"

PYTHONPATH="${BUILD_DIR}/python" python3 - <<'PY'
import vcfx

tested = {
    'alignment_checker','allele_counter','variant_counter','allele_freq_calc',
    'info_aggregator','info_parser','info_summarizer','fasta_converter',
    'allele_balance_calc','concordance_checker','genotype_query',
    'duplicate_remover','af_subsetter','missing_detector','hwe_tester',
    'inbreeding_calculator','variant_classifier','cross_sample_concordance',
    'field_extractor','ancestry_assigner','dosage_calculator',
    'ancestry_inferrer','distance_calculator','indel_normalizer','validator'
}
all_tools = vcfx.available_tools()
untested = sorted(set(all_tools) - tested)
for tool in untested:
    proc = vcfx.run_tool(tool, '--help', capture_output=True, text=True)
    assert proc.returncode == 0
print('Extra wrapper help checks OK')
PY
