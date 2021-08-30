#!/usr/bin/env python

"""
Superclasses of ChromosomeBase for generating Elements to represent 
a Chromosome structure. ChromosomeBase has functions to visualize
these chroms and to convert them to SLiM simulation commands.

These classes are for internal use only. They are exposed in user-
facing functions in :meth:`shadie.chromosome`.
"""

from typing import Union, List
import numpy as np

# internal imports
from shadie.base.elements import ElementType
from shadie.base.defaults import SYN, NONCDS, INTRON, EXON, NEUT
from shadie.chromosome.src.base_class import ChromosomeBase


class Chromosome(ChromosomeBase):
    """Builds the default shadie chromosome used for testing.

    This chromosome is exposed to users in a factory function at
    :meth:`shadie.chromosome.default`.
    """
    def __init__(self):
        super().__init__(genome_size=10001)
        self.data.loc[0] = (
            NONCDS.altname, 0, 2000, NONCDS.name, NONCDS, NONCDS.coding)
        self.data.loc[2001] = (
            EXON.altname, 2001, 4000, EXON.name, EXON, EXON.coding)
        self.data.loc[4001] = (
            INTRON.altname, 4001, 6000, INTRON.name, INTRON, INTRON.coding)
        self.data.loc[6001] = (
            EXON.altname, 6001, 8000, EXON.name, EXON, EXON.coding)
        self.data.loc[8001] = (
            NONCDS.altname, 8001, 10000, NONCDS.name, NONCDS, NONCDS.coding)

        mutations = []
        elements = [NONCDS, EXON, INTRON]
        for elem in elements:
            for mutation in elem.mlist:
                if mutation.name not in mutations:
                    mutations.append(mutation.name)
        self.mutations = mutations


class ChromosomeRandom(ChromosomeBase): 
    """Builds a random chromosome given defined element types.

    This chromosome builder is exposed to users in a factory function 
    at :meth:`shadie.chromosome.random`. 

    The chromosome will be a set length and composed randomly of 
    intron, exon, and non-cds genomic ElementTypes with their 
    relative weights scaled by args to the self.run() function. 
    Default ElementTypes are used if not entered by the user.

    Examples
    --------
    >>> chrom = ChromosomeRandom()
    >>> chrom.run()
    """
    def __init__(
        self, 
        genome_size: int=20000, 
        intron: Union[None, ElementType, List[ElementType]]=None,
        exon: Union[None, ElementType, List[ElementType]]=None,
        noncds: ElementType=None,
        seed: Union[int, None]=None,
        ):

        super().__init__(genome_size)
        self.rng = np.random.default_rng(seed)
        self.intron = intron if intron is not None else INTRON
        self.exon = exon if exon is not None else EXON
        self.noncds = noncds if noncds is not None else NONCDS
        self.genome_size = int(genome_size - 1)

        # combine all elements into a list
        elements = []
        for ele in (self.intron, self.exon, self.noncds):
            if isinstance(ele, list):
                elements += ele
            else:
                elements.append(ele)

        # extract all mutations from the elements
        mutations = []
        for elem in elements:
            for mutation in elem.mlist:
                if mutation.name not in mutations:
                    mutations.append(mutation.name)
        self.mutations = mutations

        # check the altnames of the elements... (NOTE: what is this for?)
        # for idx, ele in enumerate(self.exons):
        #     if ele.altname is None:
        #         ele.altname = f"exon-{idx + 1}"
        # for idx, ele in enumerate(self.introns):
        #     if ele.altname is None:
        #         ele.altname = f"intron-{idx + 1}"

    def get_noncds_span(self, scale:int=5000) -> int:
        """
        Draws the number of bases until the next element from an 
        exponential distribution. The scale is the average waiting
        time in number of bp.
        """
        return int(self.rng.exponential(scale=scale))

    def get_cds_spans(self, length_scale:int=1000, intron_scale:int=1000) -> List[int]:
        """
        Draws the number of exons in a fixed length space from a 
        Poisson distribution. The lam parameter is the average number
        of events per sampled region. A value of 0.005 means one 
        intron per 200bp.
        """
        cds_span = int(self.rng.exponential(scale=length_scale))
        n_introns = int(self.rng.poisson(lam=cds_span / intron_scale))
        if n_introns:
            while True:
                splits = self.rng.dirichlet(np.ones(n_introns * 2 - 1))
                splits = (splits * cds_span).astype(int)
                splits[-1] = cds_span - sum(splits[:-1])
                if all(i > 3 for i in splits):
                    break
        else:
            splits = [cds_span]
        return splits

    def run(self, noncds_scale=5000, cds_scale=1000, intron_scale=1000) -> None:
        """Generates a random chromosome by sampling elements.

        Waiting times are randomly sampled between CDS regions, as are
        the number of introns within CDS regions.
        """
        idx = 0
        while 1:
            # start with a non-cds span
            span = self.get_noncds_span(noncds_scale)
            self.data.loc[idx] = (
                self.noncds.altname, 
                idx + 1, 
                min(idx + 1 + span, self.genome_size), 
                self.noncds.name, self.noncds,
                self.noncds.coding,
            )
            idx += span + 1
            
            # get a cds span
            spans = self.get_cds_spans(cds_scale, intron_scale)

            # break if cds goes beyond the end of the genome.
            if idx + sum(spans) + len(spans) > self.genome_size:
                break

            # enter the cds into data
            for enum, span in enumerate(spans):
                # even numbered segments are the exons (0, 2, 4)
                if not enum % 2:
                    if isinstance(self.exon, list):
                        ele = self.rng.choice(self.exon)
                    else:
                        ele = self.exon
                else:
                    if isinstance(self.intron, list):                    
                        ele = self.rng.choice(self.intron)
                    else:
                        ele = self.intron
                self.data.loc[idx] = (
                    ele.altname, 
                    idx + 1, 
                    idx + span + 1, 
                    ele.name, ele,
                    ele.coding,
                )
                idx += span + 1
        self.data = self.data.sort_index()


class ChromosomeExplicit(ChromosomeBase):
    """
    Builds a chromosome dataframe from explicit instructions provided
    as start, stop positions mapped to ElementTypes. To add a regions
    without an assigned element...

    chromosome.explicit({
        (500, 1000): g1,
        (2000, 3000): g2,
        (3001, 5000): g1,
        (5000, 10000): None,
    })
    """
    def __init__(self, data):
        genome_size = 1+(max(i[1] for i in data.keys()))
        super().__init__(genome_size)

        # check data dict for proper structure
        assert all(isinstance(i, tuple) for i in data.keys()), (
            "keys of input data should be tuples of integers.")
        assert all(isinstance(i, ElementType) for i in data.values() if i), (
            "values of input data should be ElementType objects.")

        mutations = []
        for element in data.values():
            for mutation in element.mlist:
                if mutation.name not in mutations:
                    mutations.append(mutation.name)
        self.mutations = mutations

        # entere explicit dict into data
        for key in sorted(data, key=lambda x: x[0]):
            start, end = key
            if data[key] is not None:
                self.data.loc[start] = (
                    data[key].altname, start, end, data[key].name, data[key], data[key].coding)


if __name__ == "__main__":


    import shadie

    # define mutation types
    m0 = shadie.mtype(0.5, 'n', 2.0, 1.0)
    m1 = shadie.mtype(0.5, 'g', 3.0, 1.0)
    m2 = shadie.mtype(0.5, 'f', 0)
    
    # define elements types
    e0 = shadie.etype([m0, m1], [1, 2])
    e1 = shadie.etype([m2], [1])

    # # print(e0.mlist)

    # chrom = shadie.chromosome.random(100000)
    # #print(chrom.data.iloc[:50, :4])
    # #chrom.review("chromosome")

    # default = shadie.chromosome.default()
    # print(default.data)

    # design chromosome of elements
    # Do we want users to be able to put in a chromosome like this 
    # and have the gaps filled with neutral portions? YES.
    chrom = shadie.chromosome.explicit({
        (0, 500): shadie.NONCDS,
        (500, 1000): e1,
        (2000, 3000): e0,
        (3001, 5000): e1,
    })

    #elem = chrom.data.loc[500]["eltype"]
    #chrom.to_slim_mutation_types()
    test = chrom.mutations
    print(chrom.data.head())
    # chrom.inspect()
    # print(test)
    # chrom.to_slim_elements()