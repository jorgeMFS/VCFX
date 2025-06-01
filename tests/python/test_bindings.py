def test_python_bindings(vcfx):
    assert vcfx.trim("  hello  ") == "hello"
    assert vcfx.read_maybe_compressed(b"hello") == b"hello"


def test_read_file_maybe_compressed(vcfx, tmp_path):
    import gzip

    p = tmp_path / "hello.txt.gz"
    data = b"hello world\n"
    with gzip.open(p, "wb") as fh:
        fh.write(data)

    result = vcfx.read_file_maybe_compressed(str(p))
    assert result == data


def test_split_and_get_version(vcfx):
    from pathlib import Path
    import subprocess
    import sys

    assert vcfx.split("a,b,c", ",") == ["a", "b", "c"]

    root = Path(__file__).resolve().parents[2]
    expected = subprocess.check_output(
        [sys.executable, str(root / "scripts" / "extract_version.py")],
        text=True,
    ).strip()
    assert vcfx.get_version() == expected
    assert vcfx.__version__ == expected
