#!/usr/bin/env python

"""
`shadie` is a wrapper around SLiM3 that implements selection on 
alternating hapliod/diploid lifecycles and converts user-provided 
phylogeny into SLiM3-compatible subpopulation demography
"""

__version__ = "0.0.7"

from shadie.base.mutations import mtype
from shadie.base.mutations import MutationList as mlist
from shadie.base.elements import ElementType as etype
from shadie.base.elements import ElementList as elist
from shadie.sims.model import Model
from shadie import chromosome

# from shadie.base.build import 
# from shadie.chromosome import Chromosome

# from shadie.mutations import MutationList

# from shadie.elements import ElementList
# # from shadie.globals import *
# from shadie.buildchromosome import Build

# from shadie.demography import Demography
# from shadie.main import Shadie
# from shadie.postsim import PostSim

from shadie.utils import set_loglevel
set_loglevel("INFO")
