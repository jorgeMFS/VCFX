"""Python bindings for the VCFX toolkit."""

try:
    from ._vcfx import *  # type: ignore  # noqa: F401,F403
except ModuleNotFoundError:  # pragma: no cover - fallback for pure Python envs
    import gzip
    from pathlib import Path

    __all__ = [
        "trim",
        "split",
        "read_file_maybe_compressed",
        "read_maybe_compressed",
        "get_version",
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
        "AlleleFrequency",
        "InfoSummary",
        "AlleleBalance",
        "ConcordanceRow",
        "HWEResult",
        "InbreedingCoefficient",
        "VariantClassification",
        "CrossSampleConcordanceRow",
        "AncestryAssignment",
        "DosageRow",
        "AncestryInference",
        "DistanceRow",
        "IndexEntry",
    ]

    def trim(text: str) -> str:
        """Return *text* without leading/trailing whitespace."""
        return text.strip()

    def split(text: str, delim: str) -> list[str]:
        """Split *text* on *delim* returning a list."""
        return text.split(delim)

    def read_maybe_compressed(data: bytes) -> bytes:
        """Decompress *data* if it is gzip-compressed."""
        try:
            return gzip.decompress(data)
        except OSError:
            return data

    def read_file_maybe_compressed(path: str) -> bytes:
        """Read a possibly compressed file."""
        return read_maybe_compressed(Path(path).read_bytes())

    def get_version() -> str:
        """Return the toolkit version when bindings are unavailable."""
        return "0.0.0"
else:
    __all__ = [
        "trim",
        "split",
        "read_file_maybe_compressed",
        "read_maybe_compressed",
        "get_version",
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
        "AlleleFrequency",
        "InfoSummary",
        "AlleleBalance",
        "ConcordanceRow",
        "HWEResult",
        "InbreedingCoefficient",
        "VariantClassification",
        "CrossSampleConcordanceRow",
        "AncestryAssignment",
        "DosageRow",
        "AncestryInference",
        "DistanceRow",
        "IndexEntry",
    ]

__version__ = get_version()
__all__.append("__version__")
from . import tools as _tools
from typing import Callable
import subprocess

# Re-export helper functions for convenience
available_tools = _tools.available_tools
run_tool = _tools.run_tool
alignment_checker = _tools.alignment_checker
allele_counter = _tools.allele_counter
variant_counter = _tools.variant_counter
allele_freq_calc = _tools.allele_freq_calc
ancestry_assigner = _tools.ancestry_assigner
allele_balance_calc = _tools.allele_balance_calc
dosage_calculator = _tools.dosage_calculator
concordance_checker = _tools.concordance_checker
genotype_query = _tools.genotype_query
duplicate_remover = _tools.duplicate_remover
info_aggregator = _tools.info_aggregator
info_parser = _tools.info_parser
info_summarizer = _tools.info_summarizer
fasta_converter = _tools.fasta_converter
af_subsetter = _tools.af_subsetter
allele_balance_filter = _tools.allele_balance_filter
record_filter = _tools.record_filter
missing_detector = _tools.missing_detector
hwe_tester = _tools.hwe_tester
inbreeding_calculator = _tools.inbreeding_calculator
variant_classifier = _tools.variant_classifier
cross_sample_concordance = _tools.cross_sample_concordance
field_extractor = _tools.field_extractor
ancestry_inferrer = _tools.ancestry_inferrer
annotation_extractor = _tools.annotation_extractor
compressor = _tools.compressor
custom_annotator = _tools.custom_annotator
diff_tool = _tools.diff_tool
distance_calculator = _tools.distance_calculator
file_splitter = _tools.file_splitter
format_converter = _tools.format_converter
gl_filter = _tools.gl_filter
haplotype_extractor = _tools.haplotype_extractor
haplotype_phaser = _tools.haplotype_phaser
header_parser = _tools.header_parser
impact_filter = _tools.impact_filter
indel_normalizer = _tools.indel_normalizer
indexer = _tools.indexer
ld_calculator = _tools.ld_calculator
merger = _tools.merger
metadata_summarizer = _tools.metadata_summarizer
missing_data_handler = _tools.missing_data_handler
multiallelic_splitter = _tools.multiallelic_splitter
nonref_filter = _tools.nonref_filter
outlier_detector = _tools.outlier_detector
phase_checker = _tools.phase_checker
phase_quality_filter = _tools.phase_quality_filter
phred_filter = _tools.phred_filter
population_filter = _tools.population_filter
position_subsetter = _tools.position_subsetter
probability_filter = _tools.probability_filter
quality_adjuster = _tools.quality_adjuster
ref_comparator = _tools.ref_comparator
reformatter = _tools.reformatter
region_subsampler = _tools.region_subsampler
sample_extractor = _tools.sample_extractor
sorter = _tools.sorter
subsampler = _tools.subsampler
sv_handler = _tools.sv_handler
validator = _tools.validator
from .results import (
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


def __getattr__(name: str) -> Callable[..., subprocess.CompletedProcess]:
    """Return a wrapper for a VCFX command line tool.

    Parameters
    ----------
    name : str
        Name of the tool without the ``VCFX_`` prefix.

    Returns
    -------
    Callable[..., subprocess.CompletedProcess]
        Callable that runs the requested tool.

    Raises
    ------
    AttributeError
        If *name* does not correspond to an available tool.
    """
    return getattr(_tools, name)

