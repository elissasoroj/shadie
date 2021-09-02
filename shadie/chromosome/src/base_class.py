#!/usr/bin/env python

"""
The base of Chromosome type classes used to store a dataframe of 
genomic element types and with functions for converting to SLiM
code format.
"""

from typing import Optional
import itertools
import pandas as pd
from shadie.base.defaults import NONCDS
from shadie.chromosome.src.draw import (
    draw_altair_chrom_canvas_interactive,
    draw_toyplot_chrom,
)


class ChromosomeBase:
    """Base Chromosome class. 

    This class contains functions that are inherited by the Chromosome
    classes users will initialize using factory functions, such as 
    default, random, or explicit. The attributes and functions of this
    class are shared by all superclasses, such as inspect, those
    for writing to SLiM format, and the .selected_sites property.

    Parameters
    ----------
    genome_size: int
        The size of the genome in bp; bp count starts at 0
    nucleotides: bool
        Use initializeMutationTypeNuc instead of initializeMutationType
    """
    def __init__(self, 
        genome_size: int, 
        use_nucleotides: bool=False,
        ns_sites: bool =True,
        ):

        self.genome_size = int(genome_size-1)
        self.mutations = []
        self.use_nuc = use_nucleotides
        self.ns_sites = ns_sites, #turn off S/NS sites with False
        # self.ichrom = None
        self.data = pd.DataFrame(
            columns=['name', 'start', 'end', 'eltype', 'script', 'coding'],
            data=None,
        )

    def is_coding(self, idx: int=None) -> bool:
        """Return True if a genomic region is coding (includes selection).

        If idx=None this returns True/False for the whole genome.

        Parameters
        ----------
        idx: int
            The index of a genome element from `.data`.

        Returns
        -------
        bool
        """
        if idx is None:
            return self.data.coding.any()
        return bool(self.data.loc[idx].coding)

    def to_slim_mutation_types(self) -> str:
        """Returns a string with newline separated SLIM commands to 
        initialize all MutationType objects in the chromosome data.

        Example
        -------
        >>> chrom.to_slim_mutation_types()
        'initializeMutationType("m1", 0.5, "f", 0.0);\ninitializeMut...'
        """
        if self.is_coding() == 0:
            elements = []
            elements.append(NONCDS)
        else:
             elements = self.data.script.unique()
        mut_lists = [i.mlist for i in elements]
        mutations = set(itertools.chain(*mut_lists))
        return "\n  ".join([i.to_slim(nuc=self.use_nuc) for i in mutations])

    def to_slim_element_types(self) -> str:
        """Return a string with newline separated SLIM commands to 
        initialize all ElementType objects in the chromosome data.
    
        Example
        -------
        >>> chrom.to_slim_element_types()
        'initializeGenomicElementType("g1", c(m1,m2), c(3,3), mm);\ninitia...'
        """
        if self.is_coding() == 0:
            elements = []
            elements.append(NONCDS)
        else:
            elements = sorted(self.data.script.unique(), key=lambda x: x.name)
        return "\n  ".join([i.to_slim() for i in elements])

    def to_slim_elements(self) -> str:
        """Return a SLIM command string to init genomic elements.

        This newline separated command string defines the genomic 
        elements that compose the chromosome structure.

        TODO: accommodate fully neutral simulations by skipping SLiM
        and returning a TreeSequence nicely.
    
        Example
        -------
        >>> chrom.to_slim_elements()
        'initializeGenomicElement(g3, 0, 4684);\ninitializeGenomic...'
        """
        #Note: will need to fix the formatting on this chunk**
        commands = []

        if not self.ns_sites:
            for idx in self.data.index:
                ele = self.data.loc[idx]
                commands.append("initializeGenomicElement("
                    f"{ele.eltype}, {ele.start}, {ele.end});")
            return "\n  ".join(commands)

        # the entire chrom is neutral; last bp is filled with a neutral
        #genomic element, so that SLiM doesn't complain
        if self.is_coding() == 0:
            start = int(self.genome_size-1)
            end = int(self.genome_size)
            commands.extend([
                    f"initializeGenomicElement({NONCDS.name}, {start}, {end});\n",
                ])

        # iterate over int start positions of elements
        for idx in self.data.index:
            ele = self.data.loc[idx]

            # do not write this region if neutral.
            if not self.is_coding(idx):
                pass
                # commands.append(
                    # f"initializeGenomicElement({ele.eltype}, {ele.start}, {ele.end});"
                # )

            # one of them then add the initialize command to commands.
            # SYN AND NONSYN SITES
            else:
                if 1: # if writing syn vs non-syn.
                    length = ele.end - ele.start
                    commands.extend([
                        f"types = rep({ele.eltype}, asInteger(floor({length}/3)));",
                        f"starts = {ele.start} + seqLen(integerDiv({length}, 3)) * 3;",
                        "ends = starts + 1;",
                        "initializeGenomicElement(types, starts, ends);\n",
                    ])
                else:
                    # TODO
                    pass
        return "\n  ".join(commands)

    def inspect(self, width: int=700, outfile: Optional[str]=None):
        """Return an altair interactive visualization of the chromosome.

        Parameters
        ----------
        width: int
            Width of the drawing canvas in pixels.
        outfile: Optional[str]
            An optional file path to save the plot in HTML format.

        Returns
        -------
        altair.VConcatChart
            An HTML element that will display in a jupyter notebook.
        """
        return draw_altair_chrom_canvas_interactive(self, width, outfile)

    def draw(self, width: int=700, axes: Optional['toyplot.coordinates.Cartesian']=None):
        """Return a toyplot drawing of the chromosome."""
        return draw_toyplot_chrom(self, width=width, axes=axes)
