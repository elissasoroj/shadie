#!/usr/bin/env python

"""
`shadie` is a wrapper around SLiM3 that implements selection on alternating hapliod/diploid lifecycles
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography
"""

__version__ = "0.0.4"

# from shadie.shadie import Shadie
from shadie.chromosome import Chromosome
from shadie.demography import Demography
