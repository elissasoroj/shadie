#!/usr/bin/env python

"""
Install package locally for testing.
>>> conda install [deps...] -c conda-forge
>>> pip install -e . --no-deps
"""

from setuptools import setup

# build command
setup(
    name="shadie",
    version="0.2",
    author="Elissa Sorojsrisom",
    author_email="ess2239@columbia.edu",
    license="GPLv3",
    description="SLiM3 Wrapper Program, 'Simulating Haploid-Diploid Evolution'",
    install_requires = [
        "pandas",
        "numpy",
        "pyslim",
        "tskit",
        "toyplot",
        "toytree",
        "loguru",
        "altair",
    ],
    classifiers=["Programming Language :: Python :: 3"],
)
