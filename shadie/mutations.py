#!/usr/bin/env python

"""
Allows user to create mutation types for their simulation
"""

#package imports
import numpy as np
import pandas as pd
import altair as alt
from loguru import logger
import IPython


class MutationType:
    idx = 0
    """
    Makes mutations for the simulation
    """

    def __init__(
        self,
        dominance: float,
        distribution: str,
        *params: float
        ):

    
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

        #check distribution
        DISTOPTS = ['f', 'g', 'e', 'n', 'w', 's']
        if self.dist in DISTOPTS:
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
                logger.warning("'e' and 'f' distributions take 1 param")
        else: 
            pass

        if self.dist == "g" or "n" or "w":
            if len(params) == 2:
                pass
            elif len(params) != 2:
                logger.warning("'g', 'n', and 'w' distributions take 2 params")        
        else:
            pass

        #convert single param to integer insteaf of list
        # if len(params) == 1:
        #     self.distparams = params[0]
 

    def __repr__(self):
        return f"<MutationType: {self.name}, {self.dom}, {self.dist}, {self.distparams}>"

    def inspect(self):
        print(
            f"idx: {self.name}\n"
            f"dominance coefficient: {self.dom}\n"
            f"distribution: {self.dist}\n"
            f"distribution parameters: {self.distparams}\n"
            "Distribution plot:"
        )
        if self.dist == "n":
            mean = self.distparams[0]
            stddev = self.distparams[1]
            source = pd.DataFrame({"Normal Distribution": np.random.normal(self.distparams[0], self.distparams[1], 5000)})
            base = alt.Chart(source).transform_fold(
                ["Normal Distribution"],
                as_=['Mutation', 'Fitness Effect'])
            histo = base.mark_bar(
                opacity=0.3,
            ).encode(
                alt.X('Fitness Effect:Q', bin=alt.Bin(maxbins=100)),
                alt.Y('count()', stack=None))
            rule = base.mark_rule().encode(
                x='mean(Fitness Effect):Q',
                size=alt.value(2))
        elif self.dist == "g": #distparams: (mean, shape)
            if self.dist == "g":
                mean = self.distparams[0]
                shape = self.distparams[1]
                scale = mean/shape
                stddev = np.sqrt((mean**2)/shape)
                if mean < 0:
                    source = pd.DataFrame({f"{self.name}": -(np.random.gamma(shape, -scale, 5000))})
                elif mean > 0:
                    source = pd.DataFrame({f"{self.name}": np.random.gamma(shape, scale, 5000)})
            base = alt.Chart(source).transform_fold(
                [f"{self.name}"],
                as_=[f'{self.name}', 'Fitness Effect'])
            histo = base.mark_bar(
                opacity=0.3,
            ).encode(
                alt.X('Fitness Effect:Q', bin=alt.Bin(maxbins=100)),
                alt.Y('count()', stack=None))
            rule = base.mark_rule().encode(
                x='mean(Fitness Effect):Q',
                size=alt.value(2))
        elif self.dist == "w": #distparams= (scale, shape)
            mean = "NaN"
            stddev = "Nan"
            source = pd.DataFrame({"Weibull Distribution": self.distparams[0]*(np.random.weibull(self.distparams[1], 5000))})
            base = alt.Chart(source).transform_fold(
                ["Weibull Distribution"],
                as_=['Mutation', 'Fitness Effect'])
            histo = base.mark_bar(
                opacity=0.3,
            ).encode(
                alt.X('Fitness Effect:Q', bin=alt.Bin(maxbins=100)),
                alt.Y('count()', stack=None))
            rule = base.mark_rule().encode(
                x='mean(Fitness Effect):Q',
                size=alt.value(2))
        IPython.display.display_html(histo)
        print(
            f"mean: {mean}\n"
            f"standard deviation: {stddev}")



class MutationList:

    def __init__(self, *mutationtypes):

        mutationdict = {}

        for i in mutationtypes:
            if isinstance(i, MutationType):
                distlist = [str(a) for a in i.distparams]
                script = f"'{i.name}', {i.dom}, '{i.dist}', {', '.join(distlist)}"
                mutationdict[i.name] = script
        self.mutationdict = mutationdict  #dictionary of script lines
            # else:
            #     print("please enter MutationType class objects only")

        self.mutationlist = mutationtypes

    def __repr__(self):
        mutnames = []
        for mutation in self.mutationlist:
            mutnames.append(mutation.name)
        self.mutnames = mutnames
        return f"<MutationList: {self.mutnames}>"

if __name__ == "__main__":

    # generate random chromosome
    mut1 = MutationType(0.5, "f", .01)
    mut2 = MutationType(0.5, "n", .05, .02)
    print( mut1, mut2)
    list1 = MutationList(mut1, mut2)
    print(list1)
    print(list1.mutnames)

    print(list1.mutationdict)
