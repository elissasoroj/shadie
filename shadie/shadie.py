#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import subprocess
import numpy as np
import pandas as pd
import msprime
import pyslim


from shadie.globals import EXON 

from shadie.chromosome import Chromosome
from shadie.demography import Demography

class Shadie(object):
    """
    A program.script object for writing the script for simulation in SLiM3
    """
    def __init__(
        self,
        tree=None,          # reads in 
        chromosome=None,        # reads in chromosome object from Chromosome class
        Ne = 1000,          # K
        nsamples=2,         # number of sampled haplotypes per tip in final data 
        reproduction="pter",    # defines how gametes get selected and replicate
        recomb=1e-9,        # sets rate in `initializeRecombinationRate`, also accepts map
        model = "nonWF",    # nucleotide simulation must be nonWF
        treeseq = "T",      # turns on tree sequence recording 
        ):
        """
        Builds script to run SLiM3 simulation

        Parameters:
        -----------
        tree: (str)
            Optional. A newick string or Toytree object of a species tree with edges in 
            SLiM 'generation' units

        Ne (int): default = 1000
            The effective population size. This value will be set to all edges of the tree.

        organism (str):
            Options: "pteridophyte", "bryophyte", "angiosperm".
            Defines the haploid/diploid lifecycle and how selection will act at different
            life stages. Also defines how individuals replicate and hoow gametes are generated

        recomb (float):
            The per-site per-generation recombination rate.

        """

        self.recomb = recomb
        self.Ne = Ne
        self.model = model
        self.chromosome = chromosome
        self.reproduction = reproduction

        if self.chromosome == None:
            gene = Chromosome()
            self.chromosome = gene.genome 
            self.mutationlist = gene.mutationlist
            self.elementlist = gene.elementlist
            self.chromosome = gene.genome
            self.gensize = gene.gensize

        elif isinstance(self.chromosome, Chromosome): 
            self.mutationlist = self.chromosome.mutationlist
            self.elementlist = self.chromosome.elementlist
            self.chromosome = self.chromosome.genome
            self.gensize = self.chromosome.gensize

        else:
            raise ValueError("please input valid Chromosome class object")

    def write(self, filename="shadie.slim"):
        "writes the .slim script; optional to provide filename as 'filename.slim'"
        
        self.filename = filename #add check .slim extension on filename
        outname = filename[:-5] + ".trees"
        self.outname = outname
        #else:     

        #write initialize callbacks
        script = open(self.filename, "w") #overwrites old file
        init1 = (
            "initialize() {\ninitializeSLiMModelType('"+self.model+"');\n"
            f"defineConstant('K',{self.Ne});\n"
            "initializeTreeSeq();"
            "initializeSLiMOptions(nucleotideBased=T);\n"
            f"initializeAncestralNucleotides(randomNucleotides({self.gensize}));\n")
        init2 = ""
        for key in self.mutationlist.mutationdict:
            init2 += f"initializeMutationType({self.mutationlist.mutationdict[key]});\n"
            f"{key}.convertToSubstitution = T;\n"
        
        init3 = ""
        for key in self.elementlist.elementdict:
            init3 += f"initializeGenomicElementType({self.elementlist.elementdict[key]});\n"

        init4 = ""
        for index, row in self.chromosome.iterrows():
            init4 += f"initializeGenomicElement({row['eltype']}, {int(row['start'])}, {int(row['finish'])});\n"

        initfinal = "initializeRecombinationRate("+str(self.recomb)+");\n}\n"

        #######

        #write the reproduction callback
        rep1 = (
            "\nreproduction() {\n"
            "subpop.addCrossed(individual, subpop.sampleIndividuals(1));}\n")

        gens1 = (
            "\n1 early() {\n"
            "sim.addSubpop('p1', 10);}\n"
            "early() {\n"
            "p1.fitnessScaling = K / p1.individualCount;}\n"
            "\nlate() {\n"
            "inds = p1.individuals;\n"
            "catn(sim.generation + ': ' + size(inds) + ' (' + max(inds.age) + ')');}\n"
            "2000 late() {\n"
            f"sim.treeSeqOutput('{self.outname}');"
            "sim.outputFull(ages=T);}\n"
            )

        script.write(init1 + init2 + init3 + init4 + initfinal + rep1 + gens1)
        script.close

    def run(self):
        "calls SLiM to run the simulations"
        # Run the SLiM model and load the resulting .trees
        
        subprocess.check_output(["slim", "-m", "-s", "0", self.filename])
        
    def postsim(self):
        self.ts = pyslim.load(self.outname)

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
                
        print("{}\t{}\t{}".format('ancestral', 'derived', 'count'))
        for j, a in enumerate(pyslim.NUCLEOTIDES):
            for k, b in enumerate(pyslim.NUCLEOTIDES):
                print("{}\t{}\t{}".format(a, b, M[j][k]))


    def reproduction(self):
        "defines haploid/diploid life cycle offspring generation, selection, fitness effects"

        if self.reproduction == "pter":
            #define reproduction - sporophyte produces male or hermaphroditic gametophytes
            #define sex ratio?  - this may be relevant (e.g. hermaphroodite gametophytes 
            #induce formation of male-only gametophyte around them)
            #define fitness effects - selection on gametophyte
            #**no selfing with male/female??
            pass

        elif self.reproduction == "bryo":
            #define reproduction - male (sporophyte only?), female
            #define sex ratio?  - females live longer than males, 
            #define fitness effects - selection on gametophyte occurrs for longer than selection 
                # on sporophyte
            #**no selfing with male/female??
            pass

        elif self.reproduction == "angio":
            #define reproduction - hermaphrodite spermatophytes, male/female gameotphyte
            #define sex ratio?  - can sex ratio be set for haploid? male >> female for gametophyte
            #define fitness effects - selection mainly acts on sporophyte
            pass

        else:
            #base SLiM hermaphrodite? 
            pass



if __name__ == "__main__":
    pass
