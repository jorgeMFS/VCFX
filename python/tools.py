import subprocess
import shutil
import functools
import csv
import os
from typing import Any, Callable, Sequence

# Cache for storing the list of available tools once discovered
_TOOL_CACHE: list[str] | None = None

__all__ = [
    "available_tools",
    "run_tool",
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
]


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


def alignment_checker(vcf_file: str, reference: str) -> list[dict]:
    """Run ``alignment_checker`` and parse the TSV output.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file to check.
    reference : str
        Reference FASTA file.

    Returns
    -------
    list[dict]
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    rows = list(reader)

    return rows


def allele_counter(vcf_file: str, samples: Sequence[str] | None = None) -> list[dict]:
    """Run ``allele_counter`` and return allele counts.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF input file.
    samples : Sequence[str] | None, optional
        Optional subset of sample names to process.

    Returns
    -------
    list[dict]
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


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


def allele_freq_calc(vcf_file: str) -> list[dict]:
    """Calculate allele frequencies from a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF input file.

    Returns
    -------
    list[dict]
        Rows from the frequency table as dictionaries.
    """

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "allele_freq_calc",
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def ancestry_assigner(vcf_file: str, freq_file: str) -> list[dict]:
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

    rows: list[dict] = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or ":" in line:
            continue
        parts = line.split()
        if len(parts) == 2:
            rows.append({"Sample": parts[0], "Assigned_Population": parts[1]})

    return rows


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


def info_summarizer(vcf_file: str, fields: Sequence[str]) -> list[dict]:
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


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
) -> list[dict]:
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def dosage_calculator(vcf_file: str) -> list[dict]:
    """Calculate genotype dosages for each sample."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "dosage_calculator",
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def concordance_checker(
    vcf_file: str, sample1: str, sample2: str
) -> list[dict]:
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


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


def hwe_tester(vcf_file: str) -> list[dict]:
    """Run Hardy-Weinberg equilibrium test and parse TSV output."""

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "hwe_tester",
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def inbreeding_calculator(
    vcf_file: str,
    freq_mode: str = "excludeSample",
    skip_boundary: bool = False,
) -> list[dict]:
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def variant_classifier(
    vcf_file: str, append_info: bool = False
) -> list[dict] | str:
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def cross_sample_concordance(
    vcf_file: str, samples: Sequence[str] | None = None
) -> list[dict]:
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

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


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

