#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import msprime
import pyslim

class Coal(object):
	"""
    Overlay neutral mutations using pyslim
    """
    def __init__(
        self,
        shadie,
        ):
        """
        Builds script to run SLiM3 simulation

        Parameters:
        -----------
        shadie: (Shadie class object)
            Reads in Shadie class object
        """
        self.ts = shadie.ts

	def slimcoal():
	self.tscoal = pyslim.SlimTreeSequence(msprime.mutate(self.ts, rate=1e-6, keep=True))

	print(f"The tree sequence now has {tscoal.num_mutations} mutations, "
	      f"and mean pairwise nucleotide diversity is {tscoal.diversity()}.")