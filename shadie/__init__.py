#!/usr/bin/env python

"""
`shadie` is a wrapper around SLiM3 that implements selection on alternating hapliod/diploid lifecycles
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography
"""

__version__ = "0.0.6"


from shadie.mutations import MutationType
from shadie.mutations import MutationList
from shadie.elements import ElementType
from shadie.elements import ElementList
from shadie.globals import *
from shadie.buildchromosome import Build
from shadie.chromosome import Chromosome
from shadie.demography import Demography
from shadie.main import Shadie

from .utils import set_loglevel
set_loglevel("INFO")