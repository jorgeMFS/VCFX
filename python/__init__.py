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
        "allele_balance_calc",
        "concordance_checker",
        "genotype_query",
        "duplicate_remover",
        "info_aggregator",
        "info_parser",
        "info_summarizer",
        "fasta_converter",
        "AlleleFrequency",
        "InfoSummary",
        "AlleleBalance",
        "ConcordanceRow",
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
        "allele_balance_calc",
        "concordance_checker",
        "genotype_query",
        "duplicate_remover",
        "info_aggregator",
        "info_parser",
        "info_summarizer",
        "fasta_converter",
        "AlleleFrequency",
        "InfoSummary",
        "AlleleBalance",
        "ConcordanceRow",
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
allele_balance_calc = _tools.allele_balance_calc
concordance_checker = _tools.concordance_checker
genotype_query = _tools.genotype_query
duplicate_remover = _tools.duplicate_remover
info_aggregator = _tools.info_aggregator
info_parser = _tools.info_parser
info_summarizer = _tools.info_summarizer
fasta_converter = _tools.fasta_converter
from .results import (
    AlleleFrequency,
    InfoSummary,
    AlleleBalance,
    ConcordanceRow,
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

