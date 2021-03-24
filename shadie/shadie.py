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


from globals import EXON 

from .chromosome import Chromosome
from .demography import Demography

class Shadie(object):
    """
    A program.script object for writing the script for simulation in SLiM3
    """
    def __init__(
        self,
        tree=None,          # reads in 
        genome=None,        # reads in chromosome object from Chromosome class
        Ne = 1000,          # K
        nsamples=2,         # number of sampled haplotypes per tip in final data 
        reproduction="pter",    # defines how gametes get selected and replicate
        recomb=1e-9,        # sets rate in `initializeRecombinationRate`, also accepts map
        genome_size=1e6,    # will be used to calculate chromosome end (length -1)
        model = "nonWF",    # nucleotide simulation must be nonWF
        treeseq = "T",      # turns on tree sequence recording 
        chromosome = None,  # 
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
        self.genome = genome
        self.reproduction = reproduction

        if self.chromosome == None:
            self.geneltype = pandas.DataFrame(columns = ["Type", "script"], 
                rows = ["coding", EXON] )

        elif isinstance(self.chromosome, Chromosome): 
            self.chromosome.mutrate

        else:
            raise ValueError("please input valid chromosome object")

    def write(self, filename="shadie.slim"):
        "writes the .slim script; optional to provide filename as 'filename.slim'"
        pass
        
        self.filename = filename #add check .slim extension on filename


        else:
            self.mutrate = mutrate
            self.muttype = muttype
            self.geneltype = geneltype
            self.genel = genel


        #write initialize callbacks
        script = open(filename.self, "a") #appends so that user does not accidentally overwrite old simulation
        initialize = (
            "initialize() {\ninitializeSLiMModelType("+self.model+");\n"
            f"defineConstant('K',{self.Ne});\n"
            f"initializeMutationRate({self.mutrate});\n"
            f"initializeMutationType{self.muttype};\n m1.convertToSubstitution = T;\n"
            f"initializeGenomicElementType{self.geneltype};\n"
            f"initializeGenomicElement{self.genel};\n"
            "initializeRecombinationRate("+self.recomb+");\n}"
        )
        script.write(L1)
        script.close

    def simulate(self):
        "calls SLiM to run the simulations"
        # Run the SLiM model and load the resulting .trees
        pass
        
        subprocess.check_output(["slim", "-m", "-s", "0", "./recipe_17.4.slim"])
        ts = pyslim.load("./recipe_17.4.trees")


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
