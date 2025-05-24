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
from . import tools as _tools

# Re-export helper functions for convenience
available_tools = _tools.available_tools
run_tool = _tools.run_tool


def __getattr__(name):
    """Provide access to tool wrappers as attributes."""
    return getattr(_tools, name)

