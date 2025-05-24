import pathlib
import re
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

def read_version():
    """Extract the project version from the top-level CMakeLists.txt."""
    root = pathlib.Path(__file__).resolve().parent.parent / "CMakeLists.txt"
    text = root.read_text()
    major = re.search(r"set\(VCFX_VERSION_MAJOR\s+([0-9]+)\)", text)
    minor = re.search(r"set\(VCFX_VERSION_MINOR\s+([0-9]+)\)", text)
    patch = re.search(r"set\(VCFX_VERSION_PATCH\s+([0-9]+)\)", text)
    if major and minor and patch:
        return f"{major.group(1)}.{minor.group(1)}.{patch.group(1)}"
    m = re.search(r"project\(VCFX.*VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)\b", text)
    return m.group(1) if m else "0.0.0"

class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY={extdir}',
            '-DPYTHON_BINDINGS=ON'
        ]
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        source_dir = pathlib.Path(__file__).resolve().parent.parent
        subprocess.check_call(['cmake', str(source_dir)] + cmake_args, cwd=build_temp)
        subprocess.check_call(['cmake', '--build', '.', '--target', '_vcfx'], cwd=build_temp)

setup(
    name='vcfx',
    version=read_version(),
    packages=['vcfx'],
    package_dir={'vcfx': '.'},
    ext_modules=[CMakeExtension('vcfx._vcfx')],
    cmdclass={'build_ext': CMakeBuild},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
    ],
)
