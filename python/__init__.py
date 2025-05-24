"""Python bindings for the VCFX toolkit."""

from ._vcfx import *  # noqa: F401,F403
from . import tools as _tools

# Re-export helper functions for convenience
available_tools = _tools.available_tools
run_tool = _tools.run_tool


def __getattr__(name):
    """Provide access to tool wrappers as attributes."""
    return getattr(_tools, name)

