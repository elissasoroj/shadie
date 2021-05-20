#!/usr/bin/env python

"""
Allows user to create element types for their simulation

SLIM examples:
--------------
initializeGenomicElementType("g1", m1, 1.0);
initializeGenomicElementType("g1", c(m2,m3,m4), c(2,8,0.1));
initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(2.5e-5));
"""

# from loguru import logger
from typing import List
import numpy as np
import scipy.stats as stats
import toyplot
from shadie.base.mutations import MutationList, MutationTypeBase



class ElementType:
    """
    Makes element types for the simulation. An Element contains one or
    more MutationTypes and their relative frequencies.

    Parameters:
    -----------
    MutationType or MutationList (str): 
       For random generator: number of exons per gene; default = random
       For dict generator: number of exons per gene; default = 1

    frequency list(float): 
        Chance of mutation at each bp
    """    
    idx = 0
    def __init__(
        self,
        mutations: List[MutationTypeBase],
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
            if isinstance(mutations, (list, tuple, MutationTypeBase))
            else mutations
        )

        #store coding attribute (0 = noncoding)
        test = 0
        for mut in self.mlist:
            test += mut.coding
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
        view = [
            self.altname, self.name, self.mlist.names, 
            self.freq
        ]
        return f"<ElementType: {', '.join(map(str, view))}>"


    def to_slim(self):
        """
        Returns the SLIM command to define an ElementType object.
        
        Examples: 
        --------
        initializeGenomicElementType("g1", m1, 1.0);
        initializeGenomicElementType("g1", c(m2,m3,m4), c(2,8,0.1));
        initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(2.5e-5));
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
        """
        Returns a histogram of a kde mixture of points drawn from 
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
    """
    Creates a list of mutations and the mutationdict object for Shadie
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



# class ElementList:
#     """
#     Elementlist container class for multiple ElementType objects.
#     """   
    
#     def __init__(self, *elementtypes: List[ElementType]):

#         for ele in elements:
#             mutscript = [str(a) for a in i.mutations]
#             freqscript = [str(b) for b in i.freq]
#             script = (
#                 f"'{i.name}', c({', '.join(mutscript)})," 
#                 f"c({', '.join(freqscript)}), {i.mutmatrix}"
#                 )
#             elementdict[i.name] = script
#         self.elementdict = elementdict  #dictionary of script lines

#         self.elementlist = elementtypes

#     def __repr__(self):
#         elnames = []
#         for element in self.elementlist:
#             elnames.append(element.name)
#         return f"<ElementList: {elnames}>"

#     def inspect(self):
#         "plots mutation types fitness effect distributions"
#         print(
#             '\033[1m' + "Genomic Element List" + '\033[0m' + "\n"
#             f"Element types: {self.elementlist}\n"
#             f"Mutation types: {self.mutationlist}\n"
#             )
#         for element in self.elementlist:
#             print(
#                 '\033[1m' + "Genomic Element Type" + '\033[0m' + "\n"
#                 f"name: {element.name}\n"
#                 f"alternate name: {element.altname}\n"
#                 f"mutations: {element.mutations}\n"
#                 f"frequencies: {element.freq}\n"
#                 )
#             for mutation in element.muttypes:
#                 mutation.inspect()


if __name__ == "__main__":


    import shadie

    # create mutationlist
    mlist = shadie.mlist(
        shadie.mtype(0.5, 'f', 0.1),
        shadie.mtype(0.5, 'n', 0.05, 0.02),
    )

    # create elements (a mutation list with frequencies)
    ele1 = ElementType(mlist, (1, 1), )
    ele2 = ElementType(mlist, (0.5, 1))
    print(ele1)
    print(ele2)
    print(ele2.to_slim())
    print(ele2.coding)

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
