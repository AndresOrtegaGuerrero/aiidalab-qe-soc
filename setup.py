import pathlib
from setuptools import setup, find_packages
from setuptools.command.install import install

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="aiidalab-qe-soc",
    version="0.0.1",
    description="A aiidalab-qe plugin to bands/pdos including spin-orbit-coupling(SOC).",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/AndresOrtegaGuerrero/aiidalab-qe-soc",
    author="Andres Ortega-Guerrero",
    author_email="oandresg15@gmail.com",
    classifiers=[],
    packages=find_packages(),
    entry_points={
        "aiidalab_qe.properties": [
            "soc = aiidalab_qe_soc:soc",
        ],
        "aiida.workflows": [
            "soc_app.soc = aiidalab_qe_soc.socworkchain:SOCWorkChain",
        ],
    },
    install_requires=[
    ],
    package_data={},
    python_requires=">=3.6",
)

# Run post_install.py after installation
import subprocess
subprocess.run(['python', 'pseudos_install.py'])