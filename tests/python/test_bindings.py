
def test_python_bindings(vcfx):
    assert vcfx.trim("  hello  ") == "hello"
    assert vcfx.read_maybe_compressed(b"hello") == b"hello"
