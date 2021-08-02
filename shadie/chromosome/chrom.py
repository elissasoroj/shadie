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
    intron:Union[None, ElementType, list]=None,
    exon:Union[None, ElementType, list]=None,
    noncds:ElementType=None,
    intron_scale=1000,
    cds_scale=1000,
    noncds_scale=5000,    
    seed:Union[int, None]=None,
    ):
    """
    Build a chromosome of a set length composed randomly of intron,
    exon, and noncds element type regions.

    Parameters
    ----------
    genome_size: int = 20000
        The size in bp of the genome that will be generated.
    intron: ElementType = None
        An element type to represent introns. If None the default
        intron type is used: ...
    exon: ElementType = None
        An element type to represent exons. If None the default
        exon type is used: ...
    noncds: ElementType = None
        An element type to represent noncds. If None the default
        noncds type is used: ...
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


def explicit(genome_size, data):
    """
    Returns a chromosome built from a dictionary with end positions 
    (lowest starts at 0... maybe this should change to listing start 
    and end ...)  explicit instructions provided
    as start, stop positions of ElementTypes.

    chromosome.explicit({
        (500, 1000): e1,
        (2000, 3000): e0,
        (3001, 5000): e1,
    })
    """
    return ChromosomeExplicit(genome_size, data)
