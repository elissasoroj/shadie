#!/usr/bin/env python

"""
Allows user to create mutation types for their simulation
"""

#package imports
from loguru import logger


class MutationType:
    idx = 0
    """
    Makes mutations for the simulation
    """

    def __init__(
        self,
        dominance: float,
        distribution: str,
        *params: float):

    
        MutationType.idx += 1
        self.idx = MutationType.idx
        self.name = "m"+str(self.idx)

        self.dom = dominance
        self.dist = distribution
        self.distparams = params

        """
        Creates mutation types for the simulation

        Parameters:
        -----------

        dominance (float): 
           For random generator: number of exons per gene; default = random
           For dict generator: number of exons per gene; default = 1

        distribution (str): 
            Chance of mutation at each bp


        genome_size(int): default = =1e6
            Length of chromosome
        """

        disttypes = ['f', 'g', 'e', 'n', 'w', 's']
        if self.dist in disttypes:
            pass

        else:
            raise ValueError("please input valid distribution type")
            logger.info(f"Distribution type options: \n"
                "'f' = fixed fitness effect\n"
                "'g' = gamma distribution\n"
                "'e' = exponential distribution\n"
                "'n' = normal distirbution\n"
                "'w' = Weibull distribution\n"
                "'s' = Script-based distribution")

        if self.dist == "e" or "f":
            if len(params) == 1:
                pass
            elif len(params) != 1:
                pass
                
        else: 
            pass

        if self.dist == "g" or "n" or "w":
            if len(params) == 2:
                pass
            elif len(params) != 2:
                pass
                
        else:
            pass

    def __repr__(self):
        return f"<MutationType: {self.name}, {self.dom}, {self.dist}, {self.distparams}"


class MutationList:

    def __init__(self, *mutationtypes):

        for i in mutationtypes:
            if isinstance(i, MutationType):
                pass

        mutnames = []
        self.mutationlist = mutationtypes
        self.mutnames = mutnames


    def __repr__(self):
        for mutation in self.mutationlist:
            self.mutnames.append(mutation.name)
        return f"<MutationList: {self.mutnames}"

if __name__ == "__main__":

    # generate random chromosome
    mut1 = MutationType(0.5, "f", .01)
    mut2 = MutationType(0.5, "n", .05, .02)
    print( mut1, mut2)
    list1 = MutationList(mut1, mut2)
    print(list1)
    print(list1.mutnames)
