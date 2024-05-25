#!/usr/bin/env python

"""Shadie is a Python library for writing and executing SLiM scripts
and performing downstream evolutionary analyses. Shadie includes many
functions with a focus on implementing evolutionary models with
selection acting on alternating haploid/diploid lifecycles.

Example
-------
>>> import shadie
>>> with shadie.Model() as model:
>>>     model.initalize()
>>> ...
"""

__version__ = "1.0.1"

from shadie.base.defaults import NONCDS, EXON, INTRON, NEUT, DEL, BEN
from shadie.base.mutations import mtype
from shadie.base.mutation_list import MutationList as mlist
from shadie.base.elements import ElementType as etype
from shadie.base.elements import ElementList as elist
from shadie import chromosome
from shadie import reproduction
from shadie.sims.model import Model
from shadie import postsim
from shadie.utils import set_log_level

set_log_level("INFO")
