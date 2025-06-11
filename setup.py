import sys, os
from setuptools import setup
from setuptools.command.install import install
import subprocess

def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix

    if 'conda' in sys.prefix:
        return sys.prefix

    if 'mamba' in sys.prefix:
        return sys.prefix

    return None

def compile_and_install_software():
    """Used the subprocess module to compile/install the C software."""
    venv = get_virtualenv_path()
    if venv:
        path = os.path.abspath(venv) + '/bin'
    else:
        path = '/usr/local/bin'

    subprocess.check_call('mkdir -p exec && make exec/hetmers exec/extract_kmer_pairs', shell=True)
    subprocess.check_call(f'install -c exec/hetmers exec/extract_kmer_pairs ' + path, shell=True)

class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_and_install_software()
        super().run()

setup(cmdclass={'install': CustomInstall})