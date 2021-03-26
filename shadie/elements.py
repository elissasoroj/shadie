#!/usr/bin/env python

"""
Allows user to create mutation types for their simulation
"""

#package imports
from loguru import logger
from mutations import MutationList


class ElementType:
    idx = 0
    """
    Makes mutations for the simulation
    """

    def __init__(
        self,
        mutations: list(),
        frequency: list(),
        mutationoptions: MutationList):
    
        ElementType.idx += 1
        self.idx = ElementType.idx
        self.name = "g"+str(self.idx)

        self.mutations = mutations
        self.freq = frequency
        self.mutoptions = mutationoptions

        """
        Creates mutation types for the simulation

        Parameters:
        -----------

        mutations list(str): 
           For random generator: number of exons per gene; default = random
           For dict generator: number of exons per gene; default = 1

        frequency list(float): 
            Chance of mutation at each bp
        """
        if type(self.freq) is list:
            pass
        else:
            self.freq = [self.freq]
            logger.debug(f"frequencies: {self.freq}")

        if type(self.mutations) is list:
            pass
        else:
            self.mutations = [self.mutations]
            logger.debug(f"mutations: {self.mutations}")


        if len(self.mutations) == len(self.freq):
            pass

        else:
            raise TypeError("length of mutation list must = length of frequency list")

        for mutation in self.mutations:
            logger.info(f"{mutation}")
            if mutation in self.mutoptions.mutnames:
                pass  
            else:
                #raise TypeError("mutations must be in mutationoptions (MutationList object)")
                pass

    def __repr__(self):
        return f"<ElementType: {self.name}, {self.mutations}, {self.freq}"


class ElementList:
    
    def __init__(self, *elementtypes):

        for i in elementtypes:
            if isinstance(i, ElementType):
                pass

        self.elementlist = elementtypes

    def __repr__(self):
        elnames = []
        for element in self.elementlist:
            elnames.append(element.name)
        return f"<ElementList: {elnames}"

if __name__ == "__main__":

    # generate random chromosome
    from mutation import MutationType
    mut1 = MutationType(0.5, "f", .01)
    mut2 = MutationType(0.5, "n", .05, .02)
    list1 = MutationList(mut1, mut2)


    genel1 = ElementType(("m1", "m2"), (1,1), list1)
    genel2 = ElementType("m2", 1, list1)
    print(genel1, genel2)
    elemlist = ElementList(genel1, genel2)
    print(elemlist)