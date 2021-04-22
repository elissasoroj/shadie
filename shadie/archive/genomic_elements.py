#!/usr/bin/env python

"""
Base classes for chromosome segments.

SHADIE usage:
--------------
ElementType objects are used to define the types of mutations that occur
across a region of a chromosome. To generate a full chromosome of etype
objects it is easiest to use the chromosome constructor functions, such
as the explicit() dict method below, or standard(), or random().


SHADIE example:
---------------
e1 = shadie.etype([m0, m1, m2], [2, 8, 0.1])
e2 = shadie.etype([m1, m3], [9, 1])
e3 = shadie.etype([m1], [1])

chrom = shadie.chrom.explicit({
    4684: e3,
    4708: e1,
    4830: e2,
    4958: e1,
    ...
})


SLIM example:
--------------
initializeGenomicElementType(1, c(m2, m3, m4), c(2, 8, 0.1));
initializeGenomicElementType(2, c(m1, m3), c(9, 1));
initializeGenomicElementType(3, m1, 1);

initializeGenomicElement(g3, 0, 4684);
initializeGenomicElement(g1, 4685, 4708);
initializeGenomicElement(g2, 4709, 4830);
initializeGenomicElement(g1, 4831, 4958);
...
"""

from shadie.base.defaults import Neutral, Synonymous, Deleterious, Beneficial


class GenomicElement:
    idx = 0

    def __init__(self, mutations, frequencies):
        GenomicElement.idx += 1
        self.mutations = mutations
        self.frequencies = frequencies

    def __repr__(self):

        muts = ",".join([f"m{i.idx}" for i in self.mutations])
        freqs = ",".join(map(str, self.frequencies))
        components = [
            f"g{self.idx}",
            f"c({muts})",
            f"c({freqs})",
        ]
        return "(" + ", ".join(components) + ")"


class NonCoding(GenomicElement):
    """
    Only Neutral type mutations are allowed in NonCoding regions.
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end
        mutations = [Neutral(0.5, 0.0)]
        frequencies = [0.0]
        super().__init__(mutations, frequencies)


class Intron(GenomicElement):
    """
    Introns can have Neutral or Deleterious type mutations
    """
    def __init__(self, start, end, mutations, frequencies):
        super().__init__(mutations, frequencies)
        self.start = start
        self.end = end


class Exon(GenomicElement):
    """
    Exons can have Synonymous, Deleterious, or Beneficial mutations.
    """
    def __init__(self, start, end, mutations, frequencies):
        super().__init__(mutations, frequencies)
        self.start = start
        self.end = end



if __name__ == "__main__":

    seg0 = NonCoding(
        start=0, end=100,
    )

    seg1 = Intron(
        start=0, end=100, 
        mutations=[Neutral(0.5, 0.0), Deleterious(0.5, -0.3, 0.2)], 
        frequencies=[0.9, 0.1],
    )

    seg2 = Exon(
        start=0, end=100, 
        mutations=[
            Beneficial(0.5, 0.1), 
            Deleterious(0.5, -0.3, 0.2), 
            Synonymous(0.5, 0.0),
        ], 
        frequencies=[0.8, 0.1, 0.1],
    )

    print(seg0)
    print(seg1)
    print(seg2)
