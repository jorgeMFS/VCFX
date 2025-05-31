from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"


def test_tool_wrappers(vcfx):
    from vcfx.tools import analysis as ta, filters as tf
    rows = ta.alignment_checker(DATA / "align_Y.vcf", DATA / "align_refY.fa")
    assert rows and rows[0].Discrepancy_Type == "ALT_MISMATCH"

    counts = ta.allele_counter(DATA / "allele_counter_A.vcf")
    assert counts[0].Sample == "S1"
    assert isinstance(counts[0].Ref_Count, int)

    n = ta.variant_counter(DATA / "variant_counter_normal.vcf")
    assert n == 5

    freqs = ta.allele_freq_calc(DATA / "allele_freq_calc" / "simple.vcf")
    assert abs(freqs[0].Allele_Frequency - 0.5) < 1e-6

    annotated = ta.info_aggregator(DATA / "aggregator" / "basic.vcf", ["DP"])
    assert "#AGGREGATION_SUMMARY" in annotated

    parsed = ta.info_parser(DATA / "info_parser" / "basic.vcf", ["DP"])
    assert parsed[0]["DP"] == "10"

    summary = ta.info_summarizer(DATA / "info_summarizer" / "basic.vcf", ["DP"])
    assert abs(summary[0].Mean - 20.0) < 1e-6
    assert isinstance(summary[0].Median, float)

    fasta = ta.fasta_converter(DATA / "fasta_converter" / "basic.vcf")
    assert fasta.startswith(">")

    balance = ta.allele_balance_calc(DATA / "allele_balance_calc_A.vcf")
    assert abs(balance[0].Allele_Balance - 1.0) < 1e-6

    conc = ta.concordance_checker(
        DATA / "concordance_input.vcf",
        "SAMPLE1",
        "SAMPLE2",
    )
    assert conc[0].Concordance == "Concordant"

    filtered = ta.genotype_query(DATA / "genotype_query" / "sample.vcf", "0/1")
    assert filtered.startswith("##")

    dedup = ta.duplicate_remover(DATA / "allele_balance_calc_A.vcf")
    assert dedup.startswith("##")

    subset = ta.af_subsetter(DATA / "af_subsetter_A.vcf", "0.01-0.1")
    assert subset.startswith("##")

    flagged = ta.missing_detector(DATA / "concordance_missing_data.vcf")
    assert "MISSING_GENOTYPES=1" in flagged

    hwe_rows = ta.hwe_tester(DATA / "hwe_tester" / "basic_hwe.vcf")
    assert isinstance(hwe_rows[0].HWE_pvalue, float)

    coeff = ta.inbreeding_calculator(
        DATA / "inbreeding_calculator" / "single_sample_excludeSample_false.vcf",
        freq_mode="excludeSample",
    )
    assert isinstance(coeff[0].InbreedingCoefficient, float)

    classes = ta.variant_classifier(DATA / "classifier_mixed.vcf")
    assert classes[0].Classification

    xconc = ta.cross_sample_concordance(DATA / "concordance_some_mismatch.vcf")
    assert xconc[0].Concordance_Status
    assert isinstance(xconc[0].Num_Samples, int)

    fields = ta.field_extractor(DATA / "field_extractor_input.vcf", ["CHROM", "POS"])
    assert fields[0]["CHROM"] == "chr1"

    assign = ta.ancestry_assigner(
        DATA / "ancestry_assigner" / "input.vcf",
        DATA / "ancestry_assigner" / "freq.tsv",
    )
    assert assign[0].Assigned_Population == "EUR"

    dos = ta.dosage_calculator(DATA / "dosage_calculator" / "basic.vcf")
    assert dos[0].Dosages == "0,1,2"

    inf = ta.ancestry_inferrer(
        DATA / "ancestry_inferrer" / "eur_samples.vcf",
        DATA / "ancestry_inferrer" / "population_freqs.txt",
    )
    assert inf[0].Inferred_Population == "EUR"

    dist = ta.distance_calculator(DATA / "variant_counter_normal.vcf")
    assert isinstance(dist[1].DISTANCE, int)

    index_rows = ta.indexer(DATA / "indexer" / "basic.vcf")
    assert isinstance(index_rows[0].FILE_OFFSET, int)

    norm_text = ta.indel_normalizer(DATA / "basic_indel.vcf")
    assert norm_text.startswith("##")

    val_msg = tf.validator(DATA / "variant_counter_normal.vcf")
    assert "VCF" in val_msg
