from __future__ import annotations
import os
import subprocess
from pathlib import Path
import pytest

ROOT_DIR = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def build_dir(tmp_path_factory):
    build = tmp_path_factory.mktemp("pybuild")
    subprocess.check_call([
        "cmake",
        str(ROOT_DIR),
        "-DPYTHON_BINDINGS=ON",
    ], cwd=build)
    subprocess.check_call(["cmake", "--build", ".", "--parallel"], cwd=build)
    return build


@pytest.fixture()
def vcfx(build_dir, monkeypatch):
    python_dir = build_dir / "python"
    monkeypatch.syspath_prepend(str(python_dir))

    tool_dirs = []
    src_dir = build_dir / "src"
    if src_dir.is_dir():
        for sub in src_dir.glob("VCFX_*"):
            exe = sub / sub.name
            if exe.is_file():
                tool_dirs.append(str(sub))
    if tool_dirs:
        monkeypatch.setenv(
            "PATH",
            os.pathsep.join(tool_dirs + [os.environ.get("PATH", "")]),
        )

    import vcfx as module
    return module
