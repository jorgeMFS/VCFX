from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"


def test_tool_wrappers_extra(vcfx):
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
