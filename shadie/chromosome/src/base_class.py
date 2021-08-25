#!/usr/bin/env python

"""
The base of Chromosome type classes used to store a dataframe of 
genomic element types and with functions for converting to SLiM
code format.
"""

from typing import Optional
import itertools
import pandas as pd
from shadie.chromosome.src.draw import draw_altair_chrom_canvas_interactive


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
        The size of the genome in bp
    nucleotides: bool
        Use initializeMutationTypeNuc instead of initializeMutationType
    """
    def __init__(self, genome_size: int, use_nucleotides: bool=False):
        self.genome_size = int(genome_size)
        self.mutations = []
        self.use_nuc = use_nucleotides
        # self.ichrom = None
        self.data = pd.DataFrame(
            columns=['name', 'start', 'end', 'eltype', 'script', 'coding'],
            data=None,
        )

    @property
    def is_neutral(self):
        """Return True if no sites in the genome are under selection"""
        return not self.data.coding.any()

    def to_slim_mutation_types(self):
        """Returns a string with newline separated SLIM commands to 
        initialize all MutationType objects in the chromosome data.

        Example
        -------
        >>> chrom.to_slim_mutation_types()
        'initializeMutationType("m1", 0.5, "f", 0.0);\ninitializeMut...'
        """
        elements = self.data.script.unique()
        mut_lists = [i.mlist for i in elements]
        mutations = set(itertools.chain(*mut_lists))
        return "\n  ".join([i.to_slim(nuc=self.use_nuc) for i in mutations])

    def to_slim_element_types(self):
        """Return a string with newline separated SLIM commands to 
        initialize all ElementType objects in the chromosome data.
    
        Example
        -------
        >>> chrom.to_slim_element_types()
        'initializeGenomicElementType("g1", c(m1,m2), c(3,3), mm);\ninitia...'
        """
        elements = self.data.script.unique()
        return "\n  ".join([i.to_slim() for i in elements])

    def to_slim_elements(self):
        """Returns a string with newline separated SLIM commands to 
        initialize all Element objects in the chromosome data.

        TODO: accommodate fully neutral simulations by skipping SLiM
        and returning a TreeSequence nicely.
    
        Example
        -------
        >>> chrom.to_slim_elements()
        'initializeGenomicElement(g3, 0, 4684);\ninitializeGenomic...'
        """
        #Note: will need to fix the formatting on this chunk**
        commands = []

        # iterate over int start positions of elements
        for idx in self.data.index:

            # the entire chrom is neutral, do nothing, since nothing will
            # be simulated by SLiM. We will simply generate a TreeSequence
            # with metadata to return instead.
            if self.is_neutral:
                pass

            # some regions of the chrom are under selection. If this is
            # one of them then add the initialize command to commands.
            # TODO: if using a codon model then write in triplets.
            else:
                ele = self.data.loc[idx]
                commands.append(
                    f"initializeGenomicElement({ele.eltype}, {ele.start}, {ele.end});"
                )
                # TODO: COMMENTING OUT FOR NOW while working on reproduction.
                # TODO: append as separate commands to look nicer.
                # length = ele.end - ele.start
                # commands.append(
                #     f"types = rep({ele.eltype}, asInteger(floor({length}/3))); \n"
                #     f"starts = {ele.start} + seqLen(integerDiv({length}, 3)) * 3; \n   "
                #     "ends = starts + 1; \n"
                #     "initializeGenomicElement(types, starts, ends); \n"
                # )
        return "\n  ".join(commands)

    def inspect(self, width: int=700, outfile: Optional[str]=None):
        """Return an altair interactive visualization of the chromosome.

        Parameters
        ----------
        width: int
            Width of the drawing canvas in pixels.
        interactive: bool
            Allow interactive selector in the visualization.
        outfile: Optional[str]
            A file path to write the plot in HTML format.

        Returns
        -------
        altair.VConcatChart
            An HTML element that will display in a jupyter notebook.
        """
        return draw_altair_chrom_canvas_interactive(self, width, outfile)
