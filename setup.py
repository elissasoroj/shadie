
#!/usr/bin/env python

"""
Run `pip install -e .` to install local git version.
"""

import os
import re
from setuptools import setup

# parse version from init.py
with open("shadie/__init__.py") as init:
    CUR_VERSION = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        init.read(),
        re.M,
    ).group(1)

# build command
setup(
    name="shadie",
    packages=["shadie",
            "shadie.base",
            "shadie.chromosome",
            "shadie.chromosome.src",
            "shadie.postsim",
            "shadie.postsim.src",
            "shadie.reproduction",
            "shadie.sims"],
    version=CUR_VERSION,
    author="Elissa Sorojsrisom",
    author_email="ess2239@columbia.edu",
    license="GPLv3",
    description="'Simulating Haploid-Diploid Evolution', a wrapper for SLiM and msprime.",
    python_requires='>3.10.0',
    install_requires = [
        "numpy>=1.26.2",
        "pandas>=2.1.4",
        "pyslim>=1.0.3",
        "msprime>=1.3.0",
        "tskit>=0.5.6",
        "scipy>=1.11.4",
        "toytree", #v >=3.0.0
        "toyplot>=1.0.3",
        "altair>=5.2.0",
        "loguru>=0.6.0",
    ],
    classifiers=["Programming Language :: Python :: 3"],
)
