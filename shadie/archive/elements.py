#!/usr/bin/env python

"""
Allows user to create element types for their simulation
"""

#package imports
from loguru import logger
from shadie.mutations import MutationType
from shadie.mutations import MutationList


class ElementType:
    "Makes element types for the simulation"
    idx = 0

    def __init__(
        self,
        mutationtypes,
        frequency,
        mutationrate = 1e-9,
        altname = None
        ):
    
        ElementType.idx += 1
        self.idx = ElementType.idx
        self.name = "g"+str(self.idx)

        self.muttypes = mutationtypes
        self.freq = frequency
        self.mutrate = mutationrate

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
        mutations = []

        #check frequency formatting
        if isinstance(self.freq, list):
            pass
        elif isinstance(self.freq, tuple):
            self.freq = list(self.freq)
        else:
            self.freq = [self.freq]
            logger.debug(f"frequencies: {self.freq}")

        #check mutation type formatting
        if isinstance(self.muttypes, list):
            pass
        elif isinstance(self.muttypes, tuple):
            self.muttypes = list(self.muttypes)
        else:
            self.muttypes = [self.muttypes]
            logger.debug(f"mutation types: {self.muttypes}")

        #check list lengths
        if len(self.muttypes) == len(self.freq):
            pass
        else:
            raise TypeError("length of mutation list must = length of frequency list")

        #check for MutationType class and create list of names
        for mutation in self.muttypes:
            logger.debug(f"{mutation}")
            if isinstance(mutation, MutationType):
                mutations.append(mutation.name)
            elif mutation :
                raise TypeError("mutations must be MutationType object")
        
        #return list of names
        self.mutations = mutations
        self.altname = altname


    def __repr__(self):

        return (f"<ElementType: '{self.altname}', {self.name}, " 
            f"{self.mutations}, {self.freq}")


class ElementList:
    "makes elementlist and elementdict objects for Shadie"
    
    def __init__(self, mutationlist, *elementtypes):

        self.mutationlist = mutationlist

        if isinstance(self.mutationlist, MutationList):
            for element in elementtypes:
                for mutation in element.mutations:
                    if mutation in self.mutationlist.mutnames:
                        pass
                    else:
                        logger.warning(f"mutation is not in {self.mutationlist}")
                        #print(f"{mutation} is not in {self.mutationlist}")
        else:
            #logger.warning("mutationlist must be MutationList class object")
            TypeError("mutationlist must be MutationList class object")

        for i in elementtypes:
            if isinstance(i, ElementType):
                pass
            else: 
                raise TypeError("please enter ElementType class objects only")
                #logger.warning("please enter ElementType class objects only")

        elementdict = {}

        for i in elementtypes:
            mutscript = [str(a) for a in i.mutations]
            freqscript = [str(b) for b in i.freq]
            script = (
                f"'{i.name}', c({', '.join(mutscript)})," 
                f"c({', '.join(freqscript)})"
                )
            elementdict[i.name] = script
        self.elementdict = elementdict  #dictionary of script lines

        self.elementlist = elementtypes

    def __repr__(self):
        elnames = []
        for element in self.elementlist:
            elnames.append(element.name)
        return f"<ElementList: {elnames}>"

    def inspect(self):
        "plots mutation types fitness effect distributions"
        print(
            '\033[1m' + "Genomic Element List" + '\033[0m' + "\n"
            f"Element types: {self.elementlist}\n"
            f"Mutation types: {self.mutationlist}\n"
            )
        for element in self.elementlist:
            print(
                '\033[1m' + "Genomic Element Type" + '\033[0m' + "\n"
                f"name: {element.name}\n"
                f"alternate name: {element.altname}\n"
                f"mutations: {element.mutations}\n"
                f"frequencies: {element.freq}\n"
                )
            for mutation in element.muttypes:
                mutation.inspect()


if __name__ == "__main__":

    # generate random chromosome
    mut1 = MutationType(0.5, "f", .01)
    mut2 = MutationType(0.5, "n", .05, .02)
    list1 = MutationList(mut1, mut2)


    genel1 = ElementType([mut1, mut2], (1,1))
    genel2 = ElementType(mut2, 1)
    print(genel1, genel2)
    elemlist = ElementList(list1, genel1, genel2)

    from shadie import globals
    deflist = MutationList(globals.SYN, globals.DEL, globals.BEN)
    print(deflist)
    elemlist2 = ElementList(deflist, globals.EXON)
    print(elemlist2)
