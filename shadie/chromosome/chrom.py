#!/usr/bin/env python

"""
Convenience functions for constructing chromosome objects, 
and a chromosome viewer function.
"""

from typing import Union
from shadie.base.elements import ElementType
from shadie.chromosome.build import ChromosomeRandom, Chromosome, ChromosomeExplicit


def random(        
    genome_size:int=20000, 
    intron:ElementType=None,
    exon:ElementType=None,
    noncds:ElementType=None,
    intron_scale=1000,
    cds_scale=1000,
    noncds_scale=5000,    
    seed:Union[int, None]=None,
    ):
    """
    ...
    """
    # construct pandas DataFrame of ElementTypes
    elements = ChromosomeRandom(genome_size, intron, exon, noncds, seed)
    elements.run(noncds_scale, cds_scale, intron_scale)
    return elements


def default():
    """
    Returns the default 100Kb Chromosome of Elements used for simple 
    testing in shadie, and composing introns, exons, and noncds regions.
    """
    return Chromosome()


def explicit(data):
    """
    Returns a chromosome built from a dictionary with end positions 
    (lowest starts at 0... maybe this should change to listing start 
    and end ...)  explicit instructions provided
    as start, stop positions of ElementTypes.

    chromosome.explicit({
        500: g1,
        1000: g2,
        1200: g1,
        1500: g3,
        2000: g2,
    })
    """
    return ChromosomeExplicit(data)
