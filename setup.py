
#!/usr/bin/env python

"""
Call `pip install -e .` to install package locally for testing.
"""

from setuptools import setup

# build command
setup(
    name="shadie",
    version="0.0.7",
    author="Elissa Sorojsrisom",
    author_email="ess2239@columbia.edu",
    license="GPLv3",
    description="SLiM3 Wrapper Program, 'Simulating Haploid-Diploid Evolution'",
    install_requires = ["pandas", "numpy", "toyplot", "loguru", "toytree", "altair", "pyslim", "tskit"],
    classifiers=["Programming Language :: Python :: 3"],
)
