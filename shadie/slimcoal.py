#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import msprime
import pyslim

#internal imports

class Coal:
    """
    Overlay neutral mutations using pyslim
    """
    def __init__(
        self,
        ts,
        ):

        self.ts = ts

        """
        Builds script to run SLiM3 simulation

        Parameters:
        -----------
        shadie: (Shadie class object)
            Reads in Shadie class object
        """

    def slimcoal(self):
        self.tscoal = pyslim.SlimTreeSequence(msprime.mutate(self.ts, rate=1e-6, keep=True))

        print(f"The tree sequence now has {self.tscoal.num_mutations} mutations, "
          f"and mean pairwise nucleotide diversity is {self.tscoal.diversity()}.")