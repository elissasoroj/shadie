#!/usr/bin/env python

"""
Convenience functions for constructing chromosome objects.
"""

from typing import Union, Optional
from shadie.base.elements import ElementType
from shadie.chromosome.src.classes import (
    ChromosomeRandom, Chromosome, ChromosomeExplicit
)


def default(
    use_nucleotides: bool=False,
    use_synonymous_sites_in_coding: bool=True,
    ):
    """Return the default 100Kb Chromosome (nonCDS-exon-intron-exon-nonCDS).

    This chromosome type is typically used for testing purposes.

    Parameters
    ----------
    ...
    """
    return Chromosome(
        use_nucleotides, 
        use_synonymous_sites_in_coding, 
    )


def random(
    genome_size: int=20000,
    intron: Union[None, ElementType, list]=None,
    exon: Union[None, ElementType, list]=None,
    noncds: ElementType=None,
    intron_scale: int=1000,
    cds_scale: int=1000,
    noncds_scale: int=5000,
    seed: Union[int, None]=None,
    use_nucleotides: bool=False,
    use_synonymous_sites_in_coding: bool=True,
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

    Examples
    --------
    >>> chromosome = shadie.chromosom.random(genome_size=1e6)
    """
    # construct pandas DataFrame of ElementTypes
    elements = ChromosomeRandom(
        genome_size=genome_size, 
        intron=intron, 
        exon=exon, 
        noncds=noncds, 
        seed=seed,
        use_nucleotides=use_nucleotides,
        use_synonymous_sites_in_coding=use_synonymous_sites_in_coding,
    )
    elements.run(noncds_scale, cds_scale, intron_scale)
    return elements


def explicit(
    data, 
    genome_size: Optional[int]=None,
    use_nucleotides: bool=False,
    use_synonymous_sites_in_coding: bool=True,
    ):
    """Return a chromosome built from a dictionary.

    The dict should contain tuples of (start,end) positions as keys
    and shadie Element objects as values, which can be created using
    the shadie.etype factory function. Elements should begin at
    position=0 (TODO: right?)

    Parameters
    ----------
    data: dict
        A dictionary mapping tuples of (start,end) positions to
        shadie.ElementType object. Position tuples are inclusive
        on the position on both ends, and should not overlap (TODO?)

    Examples
    --------
    >>> m0 = shadie.mtype(0.2, 'n', -0.05, 0.015)  # deleterious
    >>> m1 = shadie.mtype(0.7, 'e', 0.1)           # beneficial
    >>> e0 = shadie.etype([m0, m1], [8, 0.1])
    >>> e1 = shadie.etype([m0], 1)
    >>> chromosome = chromosome.explicit({
    >>>     (0, 1000): shadie.NONCDS,
    >>>     (1001, 3000): e0,
    >>>     (3001, 5000): e1,
    >>>     (5001, 7000): shadie.NONCDS,
    >>> })
    >>> print(chromosome.data)
    """
    return ChromosomeExplicit(
        data=data,
        use_nucleotides=use_nucleotides,
        use_synonymous_sites_in_coding=use_synonymous_sites_in_coding,
    )
