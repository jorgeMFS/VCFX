from __future__ import annotations

"""Python bindings for the VCFX toolkit."""

try:
    from ._vcfx import *  # type: ignore  # noqa: F401,F403
except ModuleNotFoundError:  # pragma: no cover - fallback for pure Python envs
    import gzip
    from pathlib import Path
    from importlib.metadata import PackageNotFoundError, version

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
        try:
            return version(__package__ or "vcfx")
        except PackageNotFoundError:
            return "0.0.0"

from . import tools as _tools
from .tools import TOOL_NAMES
from typing import Callable
import subprocess
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

__version__ = get_version()

# Re-export helper functions for convenience
available_tools = _tools.available_tools
run_tool = _tools.run_tool

# Export command wrappers if the vcfx helper is available.  When the
# wrapper cannot be located on ``PATH`` importing ``vcfx`` should still
# succeed so that helper functions like :func:`trim` remain usable.  In
# that case the tool wrappers will be resolved lazily via
# ``__getattr__``.
try:  # pragma: no cover - depends on external binaries being built
    for _name in TOOL_NAMES:
        globals()[_name] = getattr(_tools, _name)
except FileNotFoundError:
    # Tools are unavailable; they will be looked up on demand.
    pass

__all__ = [
    "trim",
    "split",
    "read_file_maybe_compressed",
    "read_maybe_compressed",
    "get_version",
    "available_tools",
    "run_tool",
    *TOOL_NAMES,
    "AlignmentDiscrepancy",
    "AlleleCount",
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
    "__version__",
]


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

