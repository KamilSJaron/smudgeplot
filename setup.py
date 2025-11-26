"""
Custom build script for smudgeplot.

Compiles C binaries and bundles them inside the package for proper wheel distribution.
"""

import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info


def get_bin_dir():
    """Get the bin directory inside the package."""
    return Path(__file__).parent / "src" / "smudgeplot" / "bin"


def compile_binaries():
    """Compile the C binaries."""
    root = Path(__file__).parent
    bin_dir = get_bin_dir()
    bin_dir.mkdir(parents=True, exist_ok=True)

    # Check for C compiler
    cc = os.environ.get("CC", "gcc")
    if shutil.which(cc) is None and shutil.which("clang") is not None:
        cc = "clang"

    if shutil.which(cc) is None:
        print(f"WARNING: No C compiler found ({cc}). Skipping binary compilation.")
        print("Pre-built binaries may be included in wheel distributions.")
        return False

    cflags = os.environ.get("CFLAGS", "-O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing")

    # Add platform-specific flags
    if platform.system() == "Darwin":
        # macOS: ensure compatibility
        cflags += " -mmacosx-version-min=10.14"

    lib_dir = root / "src" / "lib"

    targets = [
        {
            "name": "hetmers",
            "sources": ["PloidyPlot.c", "libfastk.c", "matrix.c"],
            "libs": ["-lpthread", "-lm"],
        },
        {
            "name": "extract_kmer_pairs",
            "sources": ["PloidyList.c", "libfastk.c", "matrix.c"],
            "libs": ["-lpthread", "-lm"],
        },
    ]

    for target in targets:
        output = bin_dir / target["name"]
        sources = [str(lib_dir / src) for src in target["sources"]]

        cmd = [cc] + cflags.split() + ["-o", str(output)] + sources + target["libs"]

        print(f"Compiling {target['name']}...")
        print(f"  Command: {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            # Make executable
            os.chmod(output, 0o755)
            print(f"  Successfully compiled {target['name']}")
        except subprocess.CalledProcessError as e:
            print(f"  ERROR compiling {target['name']}: {e.stderr}")
            return False

    return True


class BuildPyWithBinaries(build_py):
    """Custom build_py that compiles binaries first."""

    def run(self):
        compile_binaries()
        super().run()


class DevelopWithBinaries(develop):
    """Custom develop that compiles binaries first."""

    def run(self):
        compile_binaries()
        super().run()


class EggInfoWithBinaries(egg_info):
    """Custom egg_info that compiles binaries first."""

    def run(self):
        compile_binaries()
        super().run()


setup(
    cmdclass={
        "build_py": BuildPyWithBinaries,
        "develop": DevelopWithBinaries,
        "egg_info": EggInfoWithBinaries,
    }
)
