#!/usr/bin/env python

"""

"""

#imports
import numpy as np
import pyslim

#internal imports
from shadie.main import Shadie
from shadie.slimcoal import Coal


class PostSim:
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
        if isinstance(shadie, Shadie):
            self.shadie = shadie

        self.tsraw = pyslim.load(self.shadie.outname)

        tsout = Coal(self.tsraw)
        tsout.slimcoal()

        self.tscoal = tsout.tscoal
        self.genemap = shadie.genemap

    
    def dNdS(self):
        ranges = []
        for index, row in self.genemap.iterrows():
            ranges.append(range(row['exonstart'], row['exonstop']))
            
        c_syn = []
        c_non = []
        nc_neut = []
        for mut in self.tscoal.mutations():
            if mut.derived_state != '1':
                c_non.append(mut.id)      
            elif mut.derived_state == '1':
                site = mut.site
                position = int(mut.position)
                for r in ranges:
                    if position in r:
                        c_syn.append(mut.id)
                    elif position not in r:
                        nc_neut.append(mut.id)
        
        dnds = len(c_non)/len(c_syn)
        print(f"dN/dS for the whole genome = {dnds}")


    def pointmutations(self):
        self.ts = self.tscoal
        M = [[0 for _ in pyslim.NUCLEOTIDES] for _ in pyslim.NUCLEOTIDES]
        for mut in self.ts.mutations():
            mut_list = mut.metadata["mutation_list"]
            k = np.argmax([u["slim_time"] for u in mut_list])
            derived_nuc = mut_list[k]["nucleotide"]
            if mut.parent == -1:
                acgt = self.tsraw.reference_sequence[int(mut.position)]
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
        self.ts = self.tscoal

        slim_gen = self.ts.metadata["SLiM"]["generation"]
        
        M = np.zeros((4,4,4,4), dtype='int')
        for mut in self.ts.mutations():
            pos = self.ts.site(mut.site).position
            # skip mutations at the end of the sequence
            if pos > 0 and pos < self.ts.sequence_length - 1:
                mut_list = mut.metadata["mutation_list"]
                k = np.argmax([u["slim_time"] for u in mut_list])
                derived_nuc = mut_list[k]["nucleotide"]
                left_nuc = self.ts.nucleotide_at(mut.node, pos - 1,
                time = slim_gen - mut_list[k]["slim_time"] - 1.0)
                right_nuc = self.ts.nucleotide_at(mut.node, pos + 1,
                time = slim_gen - mut_list[k]["slim_time"] - 1.0)
                if mut.parent == -1:
                    acgt = self.ts.reference_sequence[int(mut.position)]
                parent_nuc = pyslim.NUCLEOTIDES.index(acgt)
            else:
                parent_mut = self.ts.mutation(mut.parent)
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