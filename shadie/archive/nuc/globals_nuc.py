#!/usr/bin/env python

"""
global defaults for shadie
"""

from shadie.mutations import MutationType
from shadie.elements import ElementType

#saved as MutationType and ElementType objects
NEUT = MutationType(0.5, "f", 0.0)          #neutral mutation
SYN = MutationType(0.5, "f", 0.0)           #synonymous (REMOVED for now)
DEL = MutationType(0.1, "g", -0.03, 0.2)    #deleterious
BEN = MutationType(0.8, "e", 0.04)          #beneficial

EXON = ElementType([DEL, BEN], (8,0.1))             #exon
INTRON = ElementType([DEL], (1))                    #intron
NONCOD = ElementType(NEUT, 1, mutationrate = 0)     #non-coding

if __name__ == "__main__":
    pass
