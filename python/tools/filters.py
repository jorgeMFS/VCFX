from __future__ import annotations

from typing import Sequence
import csv

from .base import run_tool

__all__ = [
    "ld_calculator",
    "merger",
    "metadata_summarizer",
    "missing_data_handler",
    "multiallelic_splitter",
    "nonref_filter",
    "outlier_detector",
    "phase_checker",
    "phase_quality_filter",
    "phred_filter",
    "population_filter",
    "position_subsetter",
    "probability_filter",
    "quality_adjuster",
    "ref_comparator",
    "reformatter",
    "region_subsampler",
    "sample_extractor",
    "sorter",
    "subsampler",
    "sv_handler",
    "validator",
]


def ld_calculator(vcf_file: str, region: str | None = None) -> str:
    """Calculate pairwise linkage disequilibrium."""

    args: list[str] = []
    if region:
        args.extend(["--region", region])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "ld_calculator",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def merger(vcf_files: Sequence[str]) -> str:
    """Merge multiple VCF files."""

    args = ["--merge", ",".join(vcf_files)]

    result = run_tool(
        "merger",
        *args,
        capture_output=True,
        text=True,
    )

    return result.stdout


def metadata_summarizer(vcf_file: str) -> str:
    """Summarize metadata from a VCF file."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "metadata_summarizer",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def missing_data_handler(
    vcf_file: str,
    fill_missing: bool = False,
    default_genotype: str = "./.",
) -> str:
    """Handle or impute missing genotype data."""

    args: list[str] = []
    if fill_missing:
        args.append("--fill-missing")
        args.extend(["--default-genotype", default_genotype])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "missing_data_handler",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def multiallelic_splitter(vcf_file: str) -> str:
    """Split multi-allelic variants into bi-allelic records."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "multiallelic_splitter",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def nonref_filter(vcf_file: str) -> str:
    """Remove homozygous reference variants."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "nonref_filter",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def outlier_detector(
    vcf_file: str,
    metric: str,
    threshold: float,
    mode: str = "variant",
) -> list[dict]:
    """Detect variant or sample outliers based on metrics."""

    args = ["--metric", metric, "--threshold", str(threshold)]
    if mode == "variant":
        args.append("--variant")
    else:
        args.append("--sample")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "outlier_detector",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def phase_checker(vcf_file: str) -> str:
    """Keep only fully phased variants."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "phase_checker",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def phase_quality_filter(vcf_file: str, condition: str) -> str:
    """Filter variants by phase quality."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "phase_quality_filter",
        "--filter-pq",
        condition,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def phred_filter(
    vcf_file: str,
    threshold: float = 30.0,
    keep_missing_qual: bool = False,
) -> str:
    """Filter variants by PHRED quality."""

    args = ["--phred-filter", str(threshold)]
    if keep_missing_qual:
        args.append("--keep-missing-qual")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "phred_filter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def population_filter(vcf_file: str, pop_tag: str, pop_map: str) -> str:
    """Subset samples by population."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "population_filter",
        "--population",
        pop_tag,
        "--pop-map",
        pop_map,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def position_subsetter(vcf_file: str, region: str) -> str:
    """Extract variants within a genomic region."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "position_subsetter",
        "--region",
        region,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def probability_filter(vcf_file: str, condition: str) -> str:
    """Filter variants by genotype probability values."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "probability_filter",
        "--filter-probability",
        condition,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def quality_adjuster(
    vcf_file: str,
    func: str,
    no_clamp: bool = False,
) -> str:
    """Transform QUAL scores using a mathematical function."""

    args = ["--adjust-qual", func]
    if no_clamp:
        args.append("--no-clamp")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "quality_adjuster",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def ref_comparator(vcf_file: str, reference: str) -> str:
    """Compare variant alleles against a reference genome."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "ref_comparator",
        "--reference",
        reference,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def reformatter(
    vcf_file: str,
    compress_info: Sequence[str] | None = None,
    compress_format: Sequence[str] | None = None,
    reorder_info: Sequence[str] | None = None,
    reorder_format: Sequence[str] | None = None,
) -> str:
    """Reformat INFO and FORMAT fields."""

    args: list[str] = []
    if compress_info:
        args.extend(["--compress-info", ",".join(compress_info)])
    if compress_format:
        args.extend(["--compress-format", ",".join(compress_format)])
    if reorder_info:
        args.extend(["--reorder-info", ",".join(reorder_info)])
    if reorder_format:
        args.extend(["--reorder-format", ",".join(reorder_format)])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "reformatter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def region_subsampler(vcf_file: str, bed_file: str) -> str:
    """Subset variants based on a BED file of regions."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "region_subsampler",
        "--region-bed",
        bed_file,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def sample_extractor(vcf_file: str, samples: Sequence[str]) -> str:
    """Extract a subset of samples from a VCF."""

    args = ["--samples", " ".join(samples)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "sample_extractor",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def sorter(vcf_file: str, natural_chr: bool = False) -> str:
    """Sort variants lexicographically or using natural chromosome order."""

    args: list[str] = []
    if natural_chr:
        args.append("--natural-chr")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "sorter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def subsampler(vcf_file: str, n: int, seed: int | None = None) -> str:
    """Randomly subsample a VCF to *n* variants."""

    args = ["--subsample", str(n)]
    if seed is not None:
        args.extend(["--seed", str(seed)])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "subsampler",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def sv_handler(
    vcf_file: str,
    sv_filter_only: bool = False,
    sv_modify: bool = False,
) -> str:
    """Filter or modify structural variant annotations."""

    args: list[str] = []
    if sv_filter_only:
        args.append("--sv-filter-only")
    if sv_modify:
        args.append("--sv-modify")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "sv_handler",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def validator(
    vcf_file: str,
    strict: bool = False,
    report_dups: bool = False,
) -> str:
    """Validate a VCF file and return any messages."""

    args: list[str] = []
    if strict:
        args.append("--strict")
    if report_dups:
        args.append("--report-dups")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "validator",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout
