#!/usr/bin/env python

"""The base of Chromosome type classes used to store a dataframe of
genomic element types and with functions for converting to SLiM
code format.
"""

from typing import Optional
# import itertools
import pandas as pd
from shadie.chromosome.src.draw import (
    draw_altair_chrom_canvas_interactive,
    draw_toyplot_chrom,
)


# internal import
from shadie.base.defaults import SYN, NONCDS, INTRON, EXON, NEUT, BEN, DEL, EMPTY


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
        The size of the genome in bp; bp count starts at 0.
    use_nucleotides: bool
        Use initializeMutationTypeNuc instead of initializeMutationType.
    use_synonymous_sites_in_coding: bool
        For every three bases in a coding regions the first two will be
        of the MutationType of the coding region, but the third will be
        of a NONCODING MutationType to emulate synonmyous codon sites.

    Attributes
    -----------
    _skip_neutral_mutations: bool
        Do not simulation neutral mutations. This should be used if you
        plan to combine forward and backward simulations.
    """
    def __init__(
        self,
        genome_size: int,
        use_nucleotides: bool = False,
        use_synonymous_sites_in_coding: bool = False,
    ):
        self.genome_size = int(genome_size-1)
        """Size of the genome in discrete base pairs."""
        self.use_nucleotides = use_nucleotides
        """write MutationTypes init with use_nucleotides=True."""
        self.use_synonymous_sites_in_coding = use_synonymous_sites_in_coding
        """encodes neutral mutations at third codon sites."""
        self._skip_neutral_mutations: bool = True
        """hidden attr set in :ref:`shadie.Model.initialize`."""
        self.data = pd.DataFrame(
            columns=['name', 'start', 'end', 'eltype', 'script', 'is_coding'],
            data=None,
        )
        """DataFrame summary of the chromosome."""

    @property
    def mutations(self):
        """Return a list of all MutationType objects in the chromosome.

        If chromosome.skip_neutral_mutations=True this will skip any
        neutral mutations in the list it returns.
        """
        mutations = []
        idxs = []
        for elem in self.elements:
            for mutation in elem.mlist:
                if mutation.idx not in idxs:
                    mutations.append(mutation)
                    idxs.append(mutation.idx)
        return list(mutations)

    @property
    def elements(self):
        """Return a list of all ElementType objects in the chromosome.

        If chromosome.skip_neutral_mutations=True this will skip any
        neutral elements (contains only neutral mut) in the list.
        """
        if self._skip_neutral_mutations:
            if not self.is_coding(None):
                return list(set(self.data.script))
            else:
                return list(set(i for i in self.data.script if i.is_coding))
        else:
            return list(set(self.data.script))

    def is_coding(self, idx: int = None) -> bool:
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
            return self.data.is_coding.any()
        return bool(self.data.loc[idx].is_coding)

    def mutation_list(self) -> str:
        """Returns a string of mutation names for haploidDominanceCoeff.

        Example
        -------
        >>> chrom.mutation_list()
        'c(m1, m2).haploidDominanceCoeff=1.0;'
        """
        # gets muts
        mutations = []
        for elem in self.elements:
            for mutation in elem.mlist:
                mutations.append(mutation.name)
        unique_muts = list(set(mutations))
        return ', '.join(unique_muts)

    def to_slim_mutation_types(self) -> str:
        """Returns a string with newline separated SLIM commands to
        initialize all MutationType objects in the chromosome data.

        Example
        -------
        >>> chrom.to_slim_mutation_types()
        'initializeMutationType("m1", 0.5, "f", 0.0);\ninitializeMut...'
        """
        # gets muts and exclude neutral if _skip_neutral_mutations
        mutations = sorted(self.mutations, key=lambda x: x.name)
        return "\n  ".join([i.to_slim(nuc=self.use_nucleotides) for i in mutations])

    def to_slim_element_types(self) -> str:
        """Return a string with newline separated SLIM commands to
        initialize all ElementType objects in the chromosome data.

        Example
        -------
        >>> chrom.to_slim_element_types()
        'initializeGenomicElementType("g1", c(m1,m2), c(3,3), mm);\ninitia...'
        """
        # gets elements and exclude neutral if _skip_neutral_mutations
        elements = sorted(self.elements, key=lambda x: x.name)
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
        commands = []

        # cannot skip neutral mutations and have entirely neutral chrom.
        if self._skip_neutral_mutations & (not self.is_coding()):
            # the entire chrom is neutral; last bp is filled with a neutral
            # genomic element, so that SLiM doesn't complain
            start = int(self.genome_size-1)
            end = int(self.genome_size)
            commands.extend(
                    f"initializeGenomicElement({NONCDS.name}, {start}, {end});\n",
                )
            # raise ValueError(
            #     "Chromosome cannot have skip_neutral_mutations=True and "
            #     "be entirely neutral (non-coding).")

        # iterate over int start positions of elements
        for idx in self.data.index:
            ele = self.data.loc[idx]

            # neutral region is written depending on toggle.
            if not self.is_coding(idx):
                if not self._skip_neutral_mutations:
                    commands.append(
                        f"initializeGenomicElement({ele.eltype}, {ele.start}, {ele.end});"
                    )

            # coding region is written
            else:
                # write block as NONSYN/NONSYN/SYN triplets (codons)
                if self.use_synonymous_sites_in_coding:
                    length = ele.end - ele.start
                    commands.extend([
                        f"types = rep({ele.eltype}, asInteger(floor({length}/3)));",
                        f"starts = {ele.start} + seqLen(integerDiv({length}, 3)) * 3;",
                        "ends = starts + 1;",
                        "initializeGenomicElement(types, starts, ends);\n",
                    ])

                # write whole block as a single genomic element.
                else:
                    commands.append(
                        f"initializeGenomicElement({ele.eltype}, {ele.start}, {ele.end});"
                    )
        return "\n  ".join(commands)

    def inspect(self, width: int = 700, outfile: Optional[str] = None):
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

    def draw(self, width: int = 700, axes: Optional['toyplot.coordinates.Cartesian'] = None):
        """Return a toyplot drawing of the chromosome."""
        return draw_toyplot_chrom(self, width=width, axes=axes)


if __name__ == "__main__":
    pass
