import re
from pathlib import Path


def extract_version(cmake_path: Path | str | None = None) -> str:
    """Return the toolkit version as defined in CMakeLists.txt."""
    if cmake_path is None:
        cmake_path = Path(__file__).resolve().parents[1] / "CMakeLists.txt"
    cmake_text = Path(cmake_path).read_text()
    major = re.search(r"set\(\s*VCFX_VERSION_MAJOR\s+(\d+)\s*\)", cmake_text)
    minor = re.search(r"set\(\s*VCFX_VERSION_MINOR\s+(\d+)\s*\)", cmake_text)
    patch = re.search(r"set\(\s*VCFX_VERSION_PATCH\s+(\d+)\s*\)", cmake_text)
    if major and minor and patch:
        return f"{major.group(1)}.{minor.group(1)}.{patch.group(1)}"
    version = re.search(r"set\(\s*VCFX_VERSION\s+\"([0-9]+\.[0-9]+\.[0-9]+)\"\s*\)", cmake_text)
    if version:
        return version.group(1)
    raise RuntimeError("Unable to parse version from CMakeLists.txt")


if __name__ == "__main__":
    print(extract_version())
