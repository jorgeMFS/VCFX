#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path
import tomllib

from extract_version import extract_version


def parse_citation(path: Path) -> str:
    text = path.read_text()
    match = re.search(r"^version:\s*\"?([^\"]+)\"?", text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"Unable to find version in {path}")
    return match.group(1).strip()


def parse_pyproject(path: Path, cmake_version: str) -> str:
    with path.open("rb") as f:
        data = tomllib.load(f)
    project = data.get("project", {})
    if "version" in project:
        return str(project["version"])  # explicit version
    if "dynamic" in project and "version" in project["dynamic"]:
        dyn = data.get("tool", {}).get("setuptools", {}).get("dynamic", {}).get("version")
        if isinstance(dyn, dict) and dyn.get("attr") == "vcfx.__version__":
            return cmake_version
    raise RuntimeError("Unable to determine version from pyproject.toml")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    citation_version = parse_citation(root / "CITATION.cff")
    cmake_version = extract_version(root / "CMakeLists.txt")
    pyproject_version = parse_pyproject(root / "python" / "pyproject.toml", cmake_version)

    if citation_version == cmake_version == pyproject_version:
        print(f"Version OK: {citation_version}")
        return 0

    print("Version mismatch detected:")
    print(f"  CITATION.cff: {citation_version}")
    print(f"  pyproject.toml: {pyproject_version}")
    print(f"  CMakeLists.txt: {cmake_version}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
