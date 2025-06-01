#!/usr/bin/env python3
"""
Script to build and publish VCFX Python package to PyPI.

This script handles:
1. Building source distribution and wheels
2. Checking the package with twine
3. Publishing to TestPyPI (optional)
4. Publishing to PyPI (with confirmation)
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_command(cmd, cwd=None, check=True):
    """Run a command and handle errors."""
    print(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd, cwd=cwd, check=check, capture_output=True, text=True
        )
        if result.stdout:
            print(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        if e.stderr:
            print(f"stderr: {e.stderr}")
        if e.stdout:
            print(f"stdout: {e.stdout}")
        sys.exit(1)


def clean_build_dirs():
    """Clean previous build artifacts."""
    python_dir = Path(__file__).parent.parent / "python"
    dirs_to_clean = ["build", "dist", "*.egg-info"]

    for pattern in dirs_to_clean:
        for path in python_dir.glob(pattern):
            if path.is_dir():
                print(f"Removing {path}")
                import shutil
                shutil.rmtree(path)
            elif path.is_file():
                print(f"Removing {path}")
                path.unlink()


def build_package():
    """Build the Python package."""
    python_dir = Path(__file__).parent.parent / "python"

    print("Building package...")
    run_command([sys.executable, "-m", "build"], cwd=python_dir)

    # Check if dist directory was created
    dist_dir = python_dir / "dist"
    if not dist_dir.exists():
        print("Error: dist directory not created")
        sys.exit(1)

    files = list(dist_dir.glob("*"))
    print(f"Built files: {[f.name for f in files]}")
    return dist_dir


def check_package(dist_dir):
    """Check package with twine."""
    print("Checking package with twine...")
    run_command([sys.executable, "-m", "twine", "check", str(dist_dir / "*")])


def publish_to_testpypi(dist_dir):
    """Publish to TestPyPI."""
    print("Publishing to TestPyPI...")
    run_command([
        sys.executable, "-m", "twine", "upload",
        "--repository", "testpypi",
        str(dist_dir / "*")
    ])


def publish_to_pypi(dist_dir):
    """Publish to PyPI."""
    print("Publishing to PyPI...")
    run_command([
        sys.executable, "-m", "twine", "upload",
        str(dist_dir / "*")
    ])


def main():
    parser = argparse.ArgumentParser(
        description="Build and publish VCFX Python package"
    )
    parser.add_argument(
        "--clean", action="store_true", help="Clean build directories first"
    )
    parser.add_argument(
        "--test", action="store_true",
        help="Publish to TestPyPI instead of PyPI"
    )
    parser.add_argument(
        "--skip-build", action="store_true",
        help="Skip building, use existing dist files"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Build and check but don't publish"
    )

    args = parser.parse_args()

    python_dir = Path(__file__).parent.parent / "python"

    if args.clean:
        clean_build_dirs()

    if not args.skip_build:
        dist_dir = build_package()
    else:
        dist_dir = python_dir / "dist"
        if not dist_dir.exists():
            print(
                "Error: dist directory doesn't exist. "
                "Run without --skip-build first."
            )
            sys.exit(1)

    check_package(dist_dir)

    if args.dry_run:
        print("Dry run complete. Package built and checked successfully.")
        return

    if args.test:
        publish_to_testpypi(dist_dir)
    else:
        # Confirm before publishing to PyPI
        response = input("Publish to PyPI? This cannot be undone. (y/N): ")
        if response.lower() == 'y':
            publish_to_pypi(dist_dir)
        else:
            print("Aborted.")


if __name__ == "__main__":
    main()
