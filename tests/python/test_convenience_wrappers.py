from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"


def test_variant_counter_wrapper(vcfx):
    from vcfx.tools import analysis as ta
    count = ta.variant_counter(DATA / "variant_counter_normal.vcf")
    assert count == 5


def test_allele_freq_calc_wrapper(vcfx):
    from vcfx.tools import analysis as ta
    freqs = ta.allele_freq_calc(DATA / "allele_freq_calc" / "simple.vcf")
    assert abs(freqs[0].Allele_Frequency - 0.5) < 1e-6
