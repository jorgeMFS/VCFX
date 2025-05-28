from __future__ import annotations

import subprocess
import shutil
import functools
import csv
import os
from typing import Any, Callable, Sequence, Type, TypeVar, get_type_hints
from dataclasses import fields as dataclass_fields

from .results import (
    AlignmentDiscrepancy,
    AlleleCount,
    AlleleFrequency,
    InfoSummary,
    AlleleBalance,
    ConcordanceRow,
    HWEResult,
    InbreedingCoefficient,
    VariantClassification,
    CrossSampleConcordanceRow,
    AncestryAssignment,
    DosageRow,
    AncestryInference,
    DistanceRow,
    IndexEntry,
)

# Cache for storing the list of available tools once discovered
_TOOL_CACHE: list[str] | None = None

# List of VCFX command line tools with convenience wrappers
TOOL_NAMES: list[str] = [
    "alignment_checker",
    "allele_counter",
    "variant_counter",
    "allele_freq_calc",
    "ancestry_assigner",
    "allele_balance_calc",
    "dosage_calculator",
    "concordance_checker",
    "genotype_query",
    "duplicate_remover",
    "info_aggregator",
    "info_parser",
    "info_summarizer",
    "fasta_converter",
    "af_subsetter",
    "allele_balance_filter",
    "record_filter",
    "missing_detector",
    "hwe_tester",
    "inbreeding_calculator",
    "variant_classifier",
    "cross_sample_concordance",
    "field_extractor",
    "ancestry_inferrer",
    "annotation_extractor",
    "compressor",
    "custom_annotator",
    "diff_tool",
    "distance_calculator",
    "file_splitter",
    "format_converter",
    "gl_filter",
    "haplotype_extractor",
    "haplotype_phaser",
    "header_parser",
    "impact_filter",
    "indel_normalizer",
    "indexer",
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

# Exported names from this module
__all__ = ["available_tools", "run_tool", *TOOL_NAMES]


def available_tools(refresh: bool = False) -> list[str]:
    """Return the list of available VCFX command line tools.

    Parameters
    ----------
    refresh : bool, optional
        If ``True`` ignore any cached value and re-run ``vcfx --list``.
        Defaults to ``False``.

    Returns
    -------
    list[str]
        Names of tools discovered on ``PATH``.

    Raises
    ------
    FileNotFoundError
        If the ``vcfx`` executable cannot be found.
    """
    global _TOOL_CACHE
    if _TOOL_CACHE is not None and not refresh:
        return _TOOL_CACHE

    exe = shutil.which("vcfx")
    if exe is None:
        tools: set[str] = set()
        for path in os.environ.get("PATH", "").split(os.pathsep):
            if not path:
                continue
            try:
                for entry in os.listdir(path):
                    if entry.startswith("VCFX_"):
                        full = os.path.join(path, entry)
                        if os.path.isfile(full) and os.access(full, os.X_OK):
                            tools.add(entry[5:])
            except OSError:
                continue
        if not tools:
            raise FileNotFoundError("vcfx wrapper not found in PATH")
        _TOOL_CACHE = sorted(tools)
        return _TOOL_CACHE

    result = subprocess.run([exe, "--list"], capture_output=True, text=True)
    if result.returncode != 0:
        _TOOL_CACHE = []
    else:
        _TOOL_CACHE = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return _TOOL_CACHE


def run_tool(
    tool: str,
    *args: str,
    check: bool = True,
    capture_output: bool = False,
    text: bool = True,
    **kwargs: Any,
) -> subprocess.CompletedProcess:
    """Run a VCFX tool using :func:`subprocess.run`.

    Parameters
    ----------
    tool : str
        Name of the tool without the ``VCFX_`` prefix.
    *args : str
        Command line arguments passed to the tool.
    check : bool, optional
        If ``True`` (default) raise ``CalledProcessError`` on a non-zero
        exit status.
    capture_output : bool, optional
        Capture standard output and error and attach them to the returned
        ``CompletedProcess``. Defaults to ``False``.
    text : bool, optional
        If ``True`` decode output as text. Defaults to ``True``.
    **kwargs : Any
        Additional keyword arguments forwarded to :func:`subprocess.run`.

    Returns
    -------
    subprocess.CompletedProcess
        The completed process instance for the invoked command.

    Raises
    ------
    FileNotFoundError
        If the requested tool cannot be found on ``PATH``.
    subprocess.CalledProcessError
        If ``check`` is ``True`` and the process exits with a non-zero
        status.
    """
    exe = shutil.which(f"VCFX_{tool}")
    if exe is None:
        raise FileNotFoundError(f"VCFX tool '{tool}' not found in PATH")
    cmd = [exe, *map(str, args)]
    return subprocess.run(cmd, check=check, capture_output=capture_output, text=text, **kwargs)


def _convert_fields(rows: list[dict], converters: dict[str, Callable[[str], Any]]) -> list[dict]:
    """Apply *converters* to fields in each row in *rows*."""
    for row in rows:
        for key, func in converters.items():
            if key in row:
                try:
                    row[key] = func(row[key])
                except ValueError:
                    row[key] = float("nan")
    return rows


T = TypeVar("T")


def _tsv_to_dataclasses(
    text: str,
    cls: Type[T],
    converters: dict[str, Callable[[str], Any]] | None = None,
    fieldnames: Sequence[str] | None = None,
) -> list[T]:
    """Parse TSV *text* into instances of *cls*.

    If *converters* is not provided, the field types defined on ``cls``
    are used to cast values (``int`` and ``float``). Custom converters may
    be supplied to override this behaviour.

    Parameters
    ----------
    text : str
        TSV formatted text to parse.
    cls : Type[T]
        Dataclass type to instantiate for each row.
    converters : dict[str, Callable[[str], Any]] | None, optional
        Optional mapping of field names to converter functions.
    fieldnames : Sequence[str] | None, optional
        Explicit field names when *text* has no header row.
    """

    lines = [ln for ln in text.splitlines() if ln.strip()]
    reader = csv.DictReader(lines, delimiter="\t", fieldnames=fieldnames)
    rows = list(reader)

    if converters is None:
        converters = {}
        hints = get_type_hints(cls)
        for f in dataclass_fields(cls):  # type: ignore[arg-type]
            ftype = hints.get(f.name, f.type)
            if ftype is int:
                converters[f.name] = int
            elif ftype is float:
                converters[f.name] = float

    if converters:
        _convert_fields(rows, converters)

    return [cls(**row) for row in rows]


def alignment_checker(vcf_file: str, reference: str) -> list[AlignmentDiscrepancy]:
    """Run ``alignment_checker`` and parse the TSV output.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file to check.
    reference : str
        Reference FASTA file.

    Returns
    -------
    list[AlignmentDiscrepancy]
        Parsed rows from the discrepancy report.
    """

    result = run_tool(
        "alignment_checker",
        "--alignment-discrepancy",
        vcf_file,
        reference,
        capture_output=True,
        text=True,
    )

    return _tsv_to_dataclasses(result.stdout, AlignmentDiscrepancy)


def allele_counter(vcf_file: str, samples: Sequence[str] | None = None) -> list[AlleleCount]:
    """Run ``allele_counter`` and return allele counts.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF input file.
    samples : Sequence[str] | None, optional
        Optional subset of sample names to process.

    Returns
    -------
    list[AlleleCount]
        Parsed allele counts.
    """

    args: list[str] = []
    if samples:
        args.extend(["--samples", " ".join(samples)])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "allele_counter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, AlleleCount)


def variant_counter(vcf_file: str, strict: bool = False) -> int:
    """Count variants in a VCF using ``variant_counter``.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.
    strict : bool, optional
        If ``True`` fail on malformed lines. Defaults to ``False``.

    Returns
    -------
    int
        Number of valid variants reported by the tool.
    """

    args: list[str] = []
    if strict:
        args.append("--strict")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "variant_counter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    out = result.stdout.strip()
    if ":" in out:
        out = out.split(":", 1)[1]
    return int(out.strip())


def allele_freq_calc(vcf_file: str) -> list[AlleleFrequency]:
    """Calculate allele frequencies from a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF input file.

    Returns
    -------
    list[AlleleFrequency]
        Parsed frequency table rows.
    """

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "allele_freq_calc",
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, AlleleFrequency)


def ancestry_assigner(vcf_file: str, freq_file: str) -> list[AncestryAssignment]:
    """Assign sample ancestry using a frequency reference file."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "ancestry_assigner",
        "--assign-ancestry",
        freq_file,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(
        result.stdout,
        AncestryAssignment,
        fieldnames=["Sample", "Assigned_Population"],
    )


def info_aggregator(vcf_file: str, fields: Sequence[str]) -> str:
    """Aggregate INFO fields and return the annotated VCF text."""

    args = ["--aggregate-info", ",".join(fields)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "info_aggregator",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def info_parser(vcf_file: str, fields: Sequence[str]) -> list[dict]:
    """Parse INFO fields from a VCF file."""

    args = ["--info", ",".join(fields)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "info_parser",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def info_summarizer(vcf_file: str, fields: Sequence[str]) -> list[InfoSummary]:
    """Summarize INFO fields from a VCF file."""

    args = ["--info", ",".join(fields)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "info_summarizer",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, InfoSummary)


def fasta_converter(vcf_file: str) -> str:
    """Convert a VCF to FASTA format and return the FASTA text."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "fasta_converter",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def allele_balance_calc(
    vcf_file: str, samples: Sequence[str] | None = None
) -> list[AlleleBalance]:
    """Calculate allele balance for samples in a VCF."""

    args: list[str] = []
    if samples:
        args.extend(["--samples", " ".join(samples)])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "allele_balance_calc",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, AlleleBalance)


def dosage_calculator(vcf_file: str) -> list[DosageRow]:
    """Calculate genotype dosages for each sample."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "dosage_calculator",
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, DosageRow)


def concordance_checker(
    vcf_file: str, sample1: str, sample2: str
) -> list[ConcordanceRow]:
    """Check genotype concordance between two samples."""

    args = ["--samples", f"{sample1} {sample2}"]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "concordance_checker",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, ConcordanceRow)


def genotype_query(
    vcf_file: str, genotype: str, strict: bool = False
) -> str:
    """Filter variants by genotype pattern and return VCF text."""

    args = ["--genotype-query", genotype]
    if strict:
        args.append("--strict")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "genotype_query",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def duplicate_remover(vcf_file: str) -> str:
    """Remove duplicate variant records from a VCF."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "duplicate_remover",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def af_subsetter(vcf_file: str, af_range: str) -> str:
    """Subset variants by allele frequency range and return VCF text."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "af_subsetter",
        "--af-filter",
        af_range,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def allele_balance_filter(vcf_file: str, threshold: float) -> str:
    """Filter variants by allele balance threshold."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "allele_balance_filter",
        "--filter-allele-balance",
        str(threshold),
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def record_filter(
    vcf_file: str, criteria: str, logic: str | None = None
) -> str:
    """Filter variant records using generic expressions."""

    args = ["--filter", criteria]
    if logic:
        args.extend(["--logic", logic])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "record_filter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def missing_detector(vcf_file: str) -> str:
    """Flag variants with missing genotypes."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "missing_detector",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def hwe_tester(vcf_file: str) -> list[HWEResult]:
    """Run Hardy-Weinberg equilibrium test and parse TSV output."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "hwe_tester",
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, HWEResult)


def inbreeding_calculator(
    vcf_file: str,
    freq_mode: str = "excludeSample",
    skip_boundary: bool = False,
) -> list[InbreedingCoefficient]:
    """Compute inbreeding coefficients from a VCF."""

    args = ["--freq-mode", freq_mode]
    if skip_boundary:
        args.append("--skip-boundary")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "inbreeding_calculator",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, InbreedingCoefficient)


def variant_classifier(
    vcf_file: str, append_info: bool = False
) -> list[VariantClassification] | str:
    """Classify variants and optionally annotate the VCF."""

    args: list[str] = []
    if append_info:
        args.append("--append-info")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "variant_classifier",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    if append_info:
        return result.stdout

    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
    header_idx = None
    for i, ln in enumerate(lines):
        if ln.split("\t")[:2] == ["CHROM", "POS"]:
            header_idx = i + 1
            break
    if header_idx is not None:
        lines = lines[header_idx:]
    fieldnames = ["CHROM", "POS", "ID", "REF", "ALT", "Classification"]
    text = "\n".join(lines)
    return _tsv_to_dataclasses(text, VariantClassification, fieldnames=fieldnames)


def cross_sample_concordance(
    vcf_file: str, samples: Sequence[str] | None = None
) -> list[CrossSampleConcordanceRow]:
    """Check genotype concordance across samples."""

    args: list[str] = []
    if samples:
        args.extend(["--samples", ",".join(samples)])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "cross_sample_concordance",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, CrossSampleConcordanceRow)


def field_extractor(vcf_file: str, fields: Sequence[str]) -> list[dict]:
    """Extract fields from a VCF and return rows as dictionaries."""

    args = ["--fields", ",".join(fields)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "field_extractor",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def ancestry_inferrer(vcf_file: str, freq_file: str) -> list[AncestryInference]:
    """Infer sample ancestry using population frequencies."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "ancestry_inferrer",
        "--frequency",
        freq_file,
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, AncestryInference, None)


def annotation_extractor(vcf_file: str, fields: Sequence[str]) -> list[dict]:
    """Extract annotation fields into a table."""

    args = ["--annotation-extract", ",".join(fields)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "annotation_extractor",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def compressor(vcf_file: str, compress: bool = True) -> bytes:
    """Compress or decompress VCF data."""

    with open(vcf_file, "rb") as fh:
        inp = fh.read()

    args = ["--compress" if compress else "--decompress"]

    result = run_tool(
        "compressor",
        *args,
        capture_output=True,
        text=False,
        input=inp,
    )

    return result.stdout


def custom_annotator(vcf_file: str, annotation_file: str) -> str:
    """Add custom annotations from *annotation_file*."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "custom_annotator",
        "--add-annotation",
        annotation_file,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def diff_tool(file1: str, file2: str) -> str:
    """Compare two VCF files and return the diff text."""

    result = run_tool(
        "diff_tool",
        "--file1",
        file1,
        "--file2",
        file2,
        capture_output=True,
        text=True,
    )

    return result.stdout


def distance_calculator(vcf_file: str) -> list[DistanceRow]:
    """Calculate distances between consecutive variants."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "distance_calculator",
        capture_output=True,
        text=True,
        input=inp,
    )

    return _tsv_to_dataclasses(result.stdout, DistanceRow)


def file_splitter(
    vcf_file: str, prefix: str = "split", output_dir: str | None = None
) -> list[str]:
    """Split a VCF into per-chromosome files and return their paths."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    cwd = output_dir or os.getcwd()

    _ = run_tool(
        "file_splitter",
        "--prefix",
        prefix,
        capture_output=True,
        text=True,
        input=inp,
        cwd=cwd,
    )

    # collect generated files
    files = [
        os.path.join(cwd, f)
        for f in sorted(os.listdir(cwd))
        if f.startswith(prefix) and f.endswith(".vcf")
    ]
    return files


def format_converter(vcf_file: str, to_format: str) -> str:
    """Convert a VCF to another format (bed or csv)."""

    if to_format not in {"bed", "csv"}:
        raise ValueError("to_format must be 'bed' or 'csv'")

    arg = "--to-bed" if to_format == "bed" else "--to-csv"

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "format_converter",
        arg,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def gl_filter(vcf_file: str, condition: str, mode: str = "all") -> str:
    """Filter VCF records by genotype likelihood conditions."""

    args = ["--filter", condition]
    if mode != "all":
        args.extend(["--mode", mode])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "gl_filter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def haplotype_extractor(
    vcf_file: str,
    block_size: int = 100000,
    check_phase_consistency: bool = False,
) -> list[dict]:
    """Extract phased haplotype blocks from a VCF."""

    args = ["--block-size", str(block_size)]
    if check_phase_consistency:
        args.append("--check-phase-consistency")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "haplotype_extractor",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def haplotype_phaser(vcf_file: str, ld_threshold: float = 0.8) -> str:
    """Group variants into haplotype blocks based on LD."""

    args = ["--ld-threshold", str(ld_threshold)]

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "haplotype_phaser",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def header_parser(vcf_file: str) -> str:
    """Extract header lines from a VCF."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "header_parser",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def impact_filter(vcf_file: str, level: str) -> str:
    """Filter variants by functional impact level."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "impact_filter",
        "--filter-impact",
        level,
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def indel_normalizer(vcf_file: str) -> str:
    """Normalize indel variants by left-aligning and splitting."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "indel_normalizer",
        capture_output=True,
        text=True,
        input=inp,
    )

    return result.stdout


def indexer(vcf_file: str, header: bool = True) -> list[IndexEntry]:
    """Create a byte offset index for a VCF.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file to index.
    header : bool, optional
        Set to ``False`` if the output lacks a header row. Defaults to ``True``.
    """

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "indexer",
        capture_output=True,
        text=True,
        input=inp,
    )

    fields = ["CHROM", "POS", "FILE_OFFSET"] if not header else None

    return _tsv_to_dataclasses(result.stdout, IndexEntry, fieldnames=fields)


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


# Lazy attribute access for tool wrappers

def __getattr__(name: str) -> Callable[..., subprocess.CompletedProcess]:
    """Return a callable wrapper for a VCFX tool.

    Parameters
    ----------
    name : str
        Name of the tool as exposed on ``PATH`` without the ``VCFX_`` prefix.

    Returns
    -------
    Callable[..., subprocess.CompletedProcess]
        A function that invokes the requested tool.

    Raises
    ------
    AttributeError
        If *name* does not correspond to an available tool.
    FileNotFoundError
        If the ``vcfx`` wrapper cannot be located on ``PATH``.
    """
    tools = available_tools()
    if name in tools:
        return functools.partial(run_tool, name)
    raise AttributeError(f"module 'vcfx' has no attribute '{name}'")

