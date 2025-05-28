# mypy: ignore-errors
import pathlib
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

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
    packages=['vcfx'],
    package_dir={'vcfx': '.'},
    package_data={'vcfx': ['py.typed']},
    ext_modules=[CMakeExtension('vcfx._vcfx')],
    cmdclass={'build_ext': CMakeBuild},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
    ],
)
