from __future__ import annotations
"""Convenience wrappers for VCFX command line tools."""
import functools  # noqa: E402
from . import analysis, filters  # noqa: E402
from .base import available_tools, run_tool  # noqa: E402

globals().update({name: getattr(analysis, name) for name in analysis.__all__})
globals().update({name: getattr(filters, name) for name in filters.__all__})

# Keep a consolidated list of tool wrapper names for the public API
TOOL_NAMES: list[str] = [
    *analysis.__all__,
    *filters.__all__,
]

__all__ = ["available_tools", "run_tool", *TOOL_NAMES]


def __getattr__(name: str):
    """Return a callable wrapper for a VCFX tool."""
    if name in TOOL_NAMES:
        return functools.partial(run_tool, name)
    raise AttributeError(f"module 'vcfx.tools' has no attribute '{name}'")
