#!/usr/bin/env python

"""Allows user to create element types for their simulation.

An Element contains one or more MutationTypes and their relative 
frequencies. It is used in a simulations to simulate the stochastic
mutation process that could include mutations with different types,
each drawing from its own probability distribution of effects.

SLIM examples
-------------
>>> initializeGenomicElementType("g1", m1, 1.0);
>>> initializeGenomicElementType("g1", c(m2,m3,m4), c(2,8,0.1));
>>> initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(2.5e-5));
"""

# from loguru import logger
from typing import List
import numpy as np
from scipy import stats
import toyplot
from shadie.base.mutations import MutationType
from shadie.base.mutation_list import MutationList


class ElementType:
    """ElementType objects represent probabilities of MutationTypes
    occurring in a genomic region.

    Parameters
    ----------
    Mutations: MutationType or MutationList
        One or more MutationType objects to be used to represent a 
        genomic element.
    frequencies: List[float]: 
        Probability of each MutationType.
    altname: str
        A str name to label this ElementType.
    """    
    idx = 0
    def __init__(
        self,
        mutations: List[MutationType],
        frequencies: List[float],
        altname=None,
        ):
    
        ElementType.idx += 1
        self.idx = ElementType.idx
        self.name = f"g{self.idx}"
        self.altname = altname

        # convert to a MutationList regardless of input type
        self.mlist = (
            MutationList(*mutations) 
            if isinstance(mutations, (list, tuple, MutationType))
            else mutations
        )

        #store coding attribute (0 = noncoding)
        test = 0
        for mut in self.mlist:
            test += mut.is_coding
        if test > 0:
            self.coding = 1
        else:
            self.coding = 0

        self.freq = frequencies

        #check frequency formatting
        if isinstance(self.freq, float):
            self.freq = [self.freq]

        # check that lengths match up
        if len(self.mlist) != len(self.freq):
            raise ValueError(
                "length of mutation list must = length of frequency list")

    def __repr__(self):
        """Custom repr for ElementType object."""
        view = [
            self.altname, self.name, self.mlist.names, 
            self.freq
        ]
        return f"<ElementType: {', '.join(map(str, view))}>"


    def to_slim(self) -> str:
        """Return the SLIM command to define an ElementType object.
        
        Example SLiM code
        -----------------
        >>> initializeGenomicElementType("g1", m1, 1.0);
        >>> initializeGenomicElementType("g1", c(m2,m3,m4), c(2,8,0.1));
        >>> initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(2.5e-5));
        """
        inner = ", ".join([
            f"'{self.name}'",
            self.mlist.names[0] if len(self.mlist.names) == 1 else 
            "c({})".format(",".join(self.mlist.names)),
            str(self.freq[0]) if len(self.freq) == 1 else 
            "c({})".format(",".join(map(str, self.freq))),
        ])
        return f"initializeGenomicElementType({inner});"

    def draw(self, **kwargs):
        """Return a histogram of a kde mixture of points drawn from 
        all mutational distributions at the selected frequencies.
        """
        canvas = toyplot.Canvas(
            kwargs.get("width", 350),
            kwargs.get("height", 250)
        )
        axes = canvas.cartesian(
            xlabel="fitness effect",
            yshow=False,
        )
        mixture = np.concatenate([
            i._dist.rvs(**i._params, size=10000) * i._neg for i in self.mlist
        ])
        weights = np.concatenate([
            np.repeat(self.freq[i], 10000) for i in range(len(self.mlist))
        ])
        kernel = stats.gaussian_kde(dataset=mixture, weights=weights)
        xpoints = np.linspace(self.mlist.min, self.mlist.max, 500)
        yvalues = kernel.pdf(xpoints)
        mark = axes.fill(
            xpoints, yvalues, 
            style = {
                "stroke": toyplot.color.Palette()[0], 
                "stroke-width": 1.5,
                "fill": toyplot.color.Palette()[0],             
                "fill-opacity": 0.33,
            },           
        )
        axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
        axes.x.ticks.show = True
        return canvas, axes, mark        


class ElementList(list):
    """A list of mutations and the mutationdict object for Shadie.
    """
    def __init__(self, *elementtypes):
        super().__init__(elementtypes)

        self.dict = {i.name: i for i in self}
        self.names = [i.name for i in self]
        self.slim_dict = {}

        # for ele in elementtypes:
        #     # build a string representation of muttype for slim
        #     mut = ...
        #     params = ", ".join(map(str, mut.distparams))
        #     srep = f"'{ele.name}', {mut.dom}, '{mut.dist}', {params}"

        #     # store mutations
        #     self.dict[ele.name] = srep
        #     self.names.append(ele.name)

    def __repr__(self):
        return f"<ElementList: {self.names}>"

    @property
    def min(self):
        "min value in the 99% CI of all MutationType distributions"
        return min([i.mlist.min for i in self])

    @property
    def max(self):
        "max value in the 99% CI of all MutationType distributions"
        return max([i.mlist.max for i in self])


if __name__ == "__main__":

    import shadie

    # create mutationlist
    mlist = shadie.mlist(
        shadie.mtype(0.5, 'f', 0.1),
        shadie.mtype(0.5, 'n', 0.05, 0.02),
    )

    print(mlist)

    # create elements (a mutation list with frequencies)
    ele1 = ElementType(mlist, (1, 1), )
    ele2 = ElementType(mlist, (0.5, 1))
    print(ele1.mlist)
    print(ele2)
    print(ele2.to_slim())
    print(ele2.coding)

    for mut in ele1.mlist:
        print(mut._expr)

    # create an ElementList ()


    # noncod = ElementType(mut1, 1, altname = "nc")
    # exon1 = ElementType([mut2, mut5], [1, 1], altname = "ex1")
    # exon2 = ElementType([mut2, mut3, mut4], [9, 1, .02], altname = "ex2")
    # intron1 = ElementType([mut2, mut5], [1, 1], altname = "int1")
    # intron2 = ElementType([mut2, mut5], [1, 1], altname = "int2")

    # mycustomlist = ElementList(mutlist, noncod, exon1, exon2, intron1, intron2)    

    # create an element list
    # elelist = ElementList(ele1, ele2)


    # gen_el1 = ElementType([mut1, mut2], (1,1))
    # gen_el2 = ElementType(mut2, 1)
    # print(genel1, genel2)
    # elemlist = ElementList(list1, genel1, genel2)

    # from shadie import globals
    # deflist = MutationList(globals.SYN, globals.DEL, globals.BEN)
    # print(deflist)
    # elemlist2 = ElementList(deflist, globals.EXON)
    # print(elemlist2)
