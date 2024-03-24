#!/usr/bin/env python

"""Run `pip install -e .` to install local git version.
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
    version=CUR_VERSION,
    author="Elissa Sorojsrisom",
    author_email="ess2239@columbia.edu",
    license="GPLv3",
    description="'Simulating Haploid-Diploid Evolution', a wrapper for SLiM and msprime.",
    install_requires=[
        "pyslim>=1.0.3",
        "toytree>=3.0.0",
        "altair>=5.2.0",
    ],
    classifiers=["Programming Language :: Python :: 3"],
)
