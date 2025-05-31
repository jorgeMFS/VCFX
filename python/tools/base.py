from __future__ import annotations

import subprocess
import shutil
import csv
import os
from typing import Any, Callable, Sequence, Type, TypeVar, get_type_hints
from dataclasses import fields as dataclass_fields

# Cache for storing the list of available tools once discovered
_TOOL_CACHE: list[str] | None = None


def available_tools(refresh: bool = False) -> list[str]:
    """Return the list of available VCFX command line tools."""
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
        _TOOL_CACHE = [
            line.strip()
            for line in result.stdout.splitlines()
            if line.strip()
        ]
    return _TOOL_CACHE


def run_tool(
    tool: str,
    *args: str,
    check: bool = True,
    capture_output: bool = False,
    text: bool = True,
    **kwargs: Any,
) -> subprocess.CompletedProcess:
    """Run a VCFX tool using :func:`subprocess.run`."""
    exe = shutil.which(f"VCFX_{tool}")
    cmd: list[str]
    if exe is None:
        vcfx_wrapper = shutil.which("vcfx")
        if vcfx_wrapper is None:
            raise FileNotFoundError(f"VCFX tool '{tool}' not found in PATH")
        cmd = [vcfx_wrapper, tool, *map(str, args)]
    else:
        cmd = [exe, *map(str, args)]
    return subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=text,
        **kwargs,
    )


def _convert_fields(
    rows: list[dict],
    converters: dict[str, Callable[[str], Any]],
) -> list[dict]:
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
    """Parse TSV *text* into instances of *cls*."""
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
