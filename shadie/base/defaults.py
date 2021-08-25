#!/usr/bin/env python

"""
Global defaults for shadie. 

These default MutationTypes are used to model neutral or non-neutral 
mutations in simulation unless the user creates additional custom 
MutationTypes.
"""

from loguru import logger
from shadie.base.mutations import mtype as MutationType
from shadie.base.elements import ElementType


# MutationTypes: selection coefficient distributions
NEUT = MutationType(0.5, "f", 0.0,)      # neutral mutation
DEL = MutationType(0.1, "g", -3.0, 1.5)  # deleterious
BEN = MutationType(0.8, "e", 0.04)       # beneficial

# EXON has 80X more deleterious (-Gamma underdominant) than it has
# beneficial (+E overdominant) mutations.
EXON = ElementType(
	[DEL, BEN], 
	(8, 0.1),  
	altname="exon",
)

# INTRON has only deleterious (-Gamma) mutations
INTRON = ElementType(
	[DEL], 
	[1],  
	altname="intron",
	)

# NONCOD has only neutral (fixed) and no mutations recorded.
NONCDS = ElementType([NEUT], [1], altname="noncds")

#SYN is used by shadie for synonymous regiions inside coding regions
SYN = ElementType([NEUT], [1], altname="syn") 


if __name__ == "__main__":
    pass
