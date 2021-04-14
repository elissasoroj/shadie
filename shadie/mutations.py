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
 

    def __repr__(self):
        return (f"<MutationType: {self.name}, {self.dom}, "
            f"{self.dist}, {self.distparams}>")

    def inspect(self):
        print(
            '\033[1m' + "Mutation Type" + '\033[0m' + "\n"
            f"idx: {self.name}\n"
            f"dominance coefficient: {self.dom}\n"
            f"distribution: {self.dist}\n"
            f"distribution parameters: {self.distparams}\n"
            "Distribution plot:"
        )
        if self.dist == "f":
            print('\033[1m' + f"NONE: fixed fitness effect = {self.distparams[0]}\n" + '\033[0m')
        else:
            if self.dist == "e": #distparams = mean (lambda) = 1/scale
                mean = self.distparams[0]
                stddev = mean**2
                draws = np.random.exponential(mean, 10000)
                dist = np.histogram(draws, 50)
                source = pd.DataFrame({'counts': dist[0], 
                    'values': dist[1][:-1], 'mean':mean}, 
                    columns=['values', 'counts', 'mean'])
            elif self.dist == "n": #distparams= (mean, stddev)
                mean = self.distparams[0]
                stddev = self.distparams[1]
                draws = np.random.normal(self.distparams[0], self.distparams[1], size=10000)
            elif self.dist == "g": #distparams: (mean, shape)
                mean = self.distparams[0]
                shape = self.distparams[1]
                scale = mean/shape
                stddev = np.sqrt((mean**2)/shape)
                if mean < 0:
                    draws = -(np.random.gamma(shape, -scale, 10000))
                elif mean > 0:
                    draws = (np.random.gamma(shape, -scale, 10000))
            elif self.dist == "w": #distparams= (scale, shape)
                shape = self.distparams[1]
                scale = self.distparams[0]
                mean = scale*((shape-1)/shape)**(1/shape)
                stddev = "Nan"
                draws = self.distparams[0]*(np.random.weibull(self.distparams[1], 10000))
            dist = np.histogram(draws, 55)
            source = pd.DataFrame({'counts': dist[0], 
                'values': dist[1][:-1], 'mean':mean}, 
                columns=['values', 'counts', 'mean'])

            base = alt.Chart(source)
            histo = base.mark_bar(
                opacity=0.4,
                color = alt.ColorName("cornflowerblue")
            ).encode(
                alt.X('values:Q', axis=alt.Axis(title='Fitness Effect')),
                alt.Y('counts:Q'),
            )
            mark = base.mark_rule(color = alt.ColorName("mediumvioletred")).encode(
                x='mean:Q',
                size=alt.value(3),
                tooltip=['mean:Q']
            )    
            IPython.display.display_html(histo + mark)
            print(
                f"mean: {mean}\n"
                f"standard deviation: {stddev}\n\n\n")

class MutationList:

    def __init__(self, *mutationtypes):

        mutationdict = {}

        for i in mutationtypes:
            if isinstance(i, MutationType):
                distlist = [str(a) for a in i.distparams]
                script = f"'{i.name}', {i.dom}, '{i.dist}', {', '.join(distlist)}"
                mutationdict[i.name] = script
        self.mutationdict = mutationdict  #dictionary of script lines
        
        self.mutationlist = mutationtypes

        mutnames = []
        for mutation in self.mutationlist:
            mutnames.append(mutation.name)
        self.mutnames = mutnames

    def __repr__(self):
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
