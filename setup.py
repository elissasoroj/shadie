#!/usr/bin/env python

"""
Install package locally for testing.
>>> conda install [deps...] -c conda-forge
>>> pip install -e . --no-deps
"""

import re
from setuptools import setup

# parse version from init.py
with open("shadie/__init__.py", 'r', encoding="utf-8") as init:
    CUR_VERSION = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        init.read(),
        re.M,
    ).group(1)

# build command
setup(
    name="shadie",
    version=CUR_VERSION,
    packages=[
        "shadie",
        "shadie.base",
        "shadie.chromosome",
        "shadie.chromosome.src",
        "shadie.postsim",
        "shadie.postsim.src",
        "shadie.reproduction",
        "shadie.sims",
    ],
    author="Elissa Sorojsrisom",
    author_email="ess2239@columbia.edu",
    license="GPLv3",
    description="SLiM3 Wrapper Program, 'Simulating Haploid-Diploid Evolution'",
    install_requires = [
        "altair",
        "numpy",
        "pandas",
        "pyslim",
        "scipy",
        "toyplot",
        "toytree",
        "tskit",
    ],
    classifiers=["Programming Language :: Python :: 3"],
)
