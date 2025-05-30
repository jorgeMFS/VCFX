import sys
from pathlib import Path
import importlib.util

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("vcfx", ROOT / "python" / "__init__.py")
vcfx = importlib.util.module_from_spec(spec)
sys.modules["vcfx"] = vcfx
spec.loader.exec_module(vcfx)  # type: ignore


def test_run_tool_fallback(tmp_path, monkeypatch):
    wrapper = tmp_path / "vcfx"
    wrapper.write_text("#!/bin/sh\necho $1\n")
    wrapper.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}:/usr/bin:/bin")

    proc = vcfx.run_tool("dummytool", capture_output=True, text=True)
    assert proc.stdout.strip() == "dummytool"
