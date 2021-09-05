#!/usr/bin/env python

"""
`shadie` is a wrapper around SLiM3 that implements selection on 
alternating hapliod/diploid lifecycles and converts user-provided 
phylogeny into SLiM3-compatible subpopulation demography
"""

__version__ = "0.1.1"

from shadie.base.defaults import NONCDS, EXON, INTRON
from shadie.base.mutations import mtype
from shadie.base.mutations import MutationList as mlist
from shadie.base.elements import ElementType as etype
from shadie.base.elements import ElementList as elist
from shadie.sims.model import Model
from shadie import chromosome
from shadie import reproduction
from shadie import postsim
from shadie.utils import set_loglevel
set_loglevel("INFO")
