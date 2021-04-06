#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import numpy as np
import pandas as pd
import msprime
import pyslim

class PostSim(object):
	"""
    Post-simulation summary statistics
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

        self.shadie = shadie



	def pointmutations(self):
	        self.ts = pyslim.load(self.shadie.outname)

	        M = [[0 for _ in pyslim.NUCLEOTIDES] for _ in pyslim.NUCLEOTIDES]
	        for mut in self.ts.mutations():
	            mut_list = mut.metadata["mutation_list"]
	            k = np.argmax([u["slim_time"] for u in mut_list])
	            derived_nuc = mut_list[k]["nucleotide"]
	            if mut.parent == -1:
	                acgt = self.ts.reference_sequence[int(mut.position)]
	                parent_nuc = pyslim.NUCLEOTIDES.index(acgt)
	            else:
	                parent_mut = self.ts.mutation(mut.parent)
	                assert(parent_mut.site == mut.site)
	                parent_nuc = parent_mut.metadata["mutation_list"][0]["nucleotide"]
	            M[parent_nuc][derived_nuc] += 1
	                
	        print("{}\t{}\t{}".format('ancestr', 'derived', 'count'))
	        for j, a in enumerate(pyslim.NUCLEOTIDES):
	            for k, b in enumerate(pyslim.NUCLEOTIDES):
	                print("{}\t{}\t{}".format(a, b, M[j][k]))

	def aminoacids (self):
		self.ts = pyslim.load(self.shadie.outname)
		slim_gen = ts.metadata["SLiM"]["generation"]
		
		M = np.zeros((4,4,4,4), dtype='int')
		for mut in ts.mutations():
			pos = ts.site(mut.site).position
			# skip mutations at the end of the sequence
			if pos > 0 and pos < ts.sequence_length - 1:
				mut_list = mut.metadata["mutation_list"]
				k = np.argmax([u["slim_time"] for u in mut_list])
				derived_nuc = mut_list[k]["nucleotide"]
				left_nuc = ts.nucleotide_at(mut.node, pos - 1,
				time = slim_gen - mut_list[k]["slim_time"] - 1.0)
				right_nuc = ts.nucleotide_at(mut.node, pos + 1,
				time = slim_gen - mut_list[k]["slim_time"] - 1.0)
				if mut.parent == -1:
				acgt = ts.reference_sequence[int(mut.position)]
				parent_nuc = pyslim.NUCLEOTIDES.index(acgt)
			else:
				parent_mut = ts.mutation(mut.parent)
				assert(parent_mut.site == mut.site)
				parent_nuc = parent_mut.metadata["mutation_list"][0]["nucleotide"]
		
			M[left_nuc, parent_nuc, right_nuc, derived_nuc] += 1
		print("{}\t{}\t{}".format('ancestral', 'derived', 'count'))
		for j0, a0 in enumerate(pyslim.NUCLEOTIDES):
			for j1, a1 in enumerate(pyslim.NUCLEOTIDES):
				for j2, a2 in enumerate(pyslim.NUCLEOTIDES):
					for k, b in enumerate(pyslim.NUCLEOTIDES):
					print("{}{}{}\t{}{}{}\t{}".format(a0, a1, a2, a0, b, a2,
						M[j0, j1, j2, k]))