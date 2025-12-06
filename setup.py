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
from setuptools.dist import Distribution

class CompilationError(Exception):
    """Raised when C binary compilation fails."""
    pass


def get_bin_dir():
    """Get the bin directory inside the package."""
    return Path(__file__).parent / "src" / "smudgeplot" / "bin"


def compile_binaries(allow_failure=False):
    """
    Compile the C binaries.

    Args:
        allow_failure: If True, missing compiler is not fatal (for editable installs).
                      Compilation errors are always fatal.

    Raises:
        CompilationError: If compilation fails or no compiler found (when allow_failure=False)
    """
    root = Path(__file__).parent
    bin_dir = get_bin_dir()
    bin_dir.mkdir(parents=True, exist_ok=True)

    # Check for C compiler
    cc = os.environ.get("CC", "gcc")
    if shutil.which(cc) is None and shutil.which("clang") is not None:
        cc = "clang"

    if shutil.which(cc) is None:
        msg = (
            f"No C compiler found (tried: {cc}, clang). "
            "Please install gcc or clang to build smudgeplot.\n"
            "On macOS: xcode-select --install\n"
            "On Ubuntu/Debian: apt-get install build-essential\n"
            "On Fedora/RHEL: dnf install gcc"
        )
        if allow_failure:
            print(f"WARNING: {msg}")
            print("Binaries will not be available until you reinstall with a compiler.")
            return False
        raise CompilationError(msg)

    cflags = os.environ.get("CFLAGS", "-O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing")

    # Add platform-specific flags
    if platform.system() == "Darwin":
        # macOS: ensure compatibility with older versions
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
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            # Make executable
            os.chmod(output, 0o755)
            print(f"  Successfully compiled {target['name']}")
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to compile {target['name']}:\n{e.stderr}"
            print(f"  ERROR: {error_msg}", file=sys.stderr)
            raise CompilationError(error_msg) from e

    # Verify all binaries exist
    for target in targets:
        binary = bin_dir / target["name"]
        if not binary.exists():
            raise CompilationError(f"Binary {target['name']} was not created")

    print(f"All binaries compiled successfully to {bin_dir}")
    return True


class BuildPyWithBinaries(build_py):
    """Custom build_py that compiles binaries first (required for wheel builds)."""

    def run(self):
        # For wheel builds, compilation must succeed
        compile_binaries(allow_failure=False)
        super().run()


class DevelopWithBinaries(develop):
    """Custom develop that compiles binaries (allows failure for editable installs)."""

    def run(self):
        # For editable installs, allow missing compiler (user can fix later)
        try:
            compile_binaries(allow_failure=True)
        except CompilationError as e:
            print(f"WARNING: {e}", file=sys.stderr)
            print("Continuing with editable install, but binaries will not work.", file=sys.stderr)
        super().run()

class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

setup(
    cmdclass={
        "build_py": BuildPyWithBinaries,
        "develop": DevelopWithBinaries,
    },
    distclass=BinaryDistribution
)
