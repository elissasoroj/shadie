#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#package imports
import io
import pandas as pd
import numpy as np
from loguru import logger

#internal imports
from shadie.mutations import MutationList
from shadie.elements import ElementList

#default element types
from shadie.globals import NONCOD
from shadie.globals import INTRON
from shadie.globals import EXON

#default mutation types
from shadie.globals import NEUT
from shadie.globals import SYN
from shadie.globals import DEL
from shadie.globals import BEN

#optional imports
try:
    import IPython 
except ImportError:
    pass

#


class Build:
    """
    Chromosome object created based on user parameters
    """

    def __init__(
    	self,
    	chromtype = "random",       #"random" or "dict" 
    	genes=None,             #number of genes on the chromosome (if None, random)
		exons = None,           #number of exons per gene (if None, random)
		mutationtypes = None,   #read in MutationList object
		elementtypes = None,    #read in ElementList object
		genome_size=1e6,        #will be used to calculate chromosome end (length -1)     
        ):
    
        self.type = chromtype
        self.genes = genes
        self.exons = exons
        self.gensize = genome_size

        """
        Builds the chromosome for SLiM3 simulation

        Parameters:
        -----------
        chromtype(str): default = "random"
            "random" will generate a random chromosome
            "dict" will accept a dictionary object

        genes(int): default = None
            Number of genes on chromosome
            Can be supplied as `int` for random chromosome generation
            or as dict object for explicit chromosome generation, e.g:

            dict(name = "name", mutations = (fitness effect,  dominance), start = bp, end = bp)
            dict(name = "a", mutations = (-.01, 0.5), start = 2000, end = 3000)

            # call simulate with details on genome structure (and which life stage selection occurs?)
            mod.simulate(
            dict(name='a', selection=-0.01, start=2000, end=3000, haploid=True),
            dict(name='b', selection=-0.02, start=5000, end=6000, haploid=True),
            dict(name='c', selection=-0.03, start=9000, end=10000, haploid=True),
            dict(name='A', selection=0.01, start=2000, end=3000, haploid=False),
            dict(name='B', selection=0.02, start=5000, end=6000, haploid=False),
            dict(name='C', selection=0.03, start=9000, end=10000, haploid=False),
            )

        introns (int): default = None
            For random generator: number of introns per gene; default = random
            For dict generator: number of introns per gene; default = 0

        exons (int): default = None
           For random generator: number of exons per gene; default = random
           For dict generator: number of exons per gene; default = 1

        mutation_rate(int): default = 1e-7
            Chance of mutation at each bp


        genome_size(int): default = =1e6
            Length of chromosome
        """

        if self.type == "custom":
            if isinstance(args, io.TextIOBase):
                self.genedf = pd.read_csv(args)

            elif isinstance(args, pd.DataFrame):
                self.genedf = args

            else:
                genelist = []
                for i in args:
                    if isinstance(i, dict):
                        genelist.append(i)
                self.genedf = pd.DataFrame(genelist)

        elif self.type == "random":
            pass
        

    def dict(self):
        "generates chromosome based on explicit user structure"
        for gene in self.genedf:
            #make the mutation types
            mutdict = {}
            #...
            self.mutdict = mutdict #my linter is telling me these should be defined in init

            #make the genomic element types
            elemdict = {}
            #...
            self.selemdict = elemdict

            #write the initializeGenomicElement lines to a dict
            initdict = {}
            #...
            self.initdict = initdict


    def random(self):
        "generates a random chromosome"
        if self.type == "random":
            if self.exons == None and self.introns == None:
                genelements = pd.DataFrame(
                columns=['name', 'start', 'finish', 'eltype', 'script'],
                data=None,
                )
                base = int(0)
                finalnc_length = np.random.randint(100, 5000)
                end = self.gensize - finalnc_length
                logger.debug("Made objects: base = {base}, genelements = {}, end = {end}")

                while base < end:
                    #make initial non-coding region
                    nc_length = np.random.randint(100, 5000)
                    genelements.loc[base, 'name'] = "noncoding"
                    genelements.loc[base, 'eltype'] = NONCOD.name
                    genelements.loc[base, 'script'] = NONCOD
                    genelements.loc[base, 'start'] = base
                    genelements.loc[base, 'finish'] = base + nc_length - 1
                    base = base + nc_length
                
                    #make first exon
                    ex_length = round(np.random.lognormal(np.log(250), np.log(1.3))) + 1
                    genelements.loc[base, 'name'] = "exon"
                    genelements.loc[base, 'eltype'] = EXON.name
                    genelements.loc[base, 'script'] = EXON
                    genelements.loc[base, 'start'] = base
                    genelements.loc[base, 'finish'] = base + ex_length -1
                    base = base + ex_length
                    logger.info("Gene added")    
                    
                    while np.random.random_sample() < 0.75:  #25% probability of stopping
                        in_length = round(np.random.normal(450, 100))
                        genelements.loc[base, 'name'] = "intron"
                        genelements.loc[base, 'eltype'] = INTRON.name
                        genelements.loc[base, 'script'] = INTRON
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + in_length -1
                        base = base + in_length
                      
                        ex_length = round(np.random.lognormal(np.log(250), np.log(1.3))) + 1
                        #do you need math for that? You can just do it with numpy
                        genelements.loc[base, 'name'] = "exon"
                        genelements.loc[base, 'eltype'] = EXON.name
                        genelements.loc[base, 'script'] = EXON
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + ex_length -1
                        base = base + ex_length 
                          
            #final non-coding region
            genelements.loc[base, 'name'] = "noncoding"
            genelements.loc[base, 'eltype'] = NONCOD.name
            genelements.loc[base, 'script'] = NONCOD
            genelements.loc[base, 'start'] = base
            genelements.loc[base, 'finish'] = self.gensize - 1
            logger.info("Chromosome complete!")
            logger.debug(genelements)
            self.genelements = genelements
            self.mutationlist = MutationList(NEUT, SYN, DEL, BEN)
            self.elementlist = ElementList(EXON, INTRON, NONCOD)

        elif self.type == "dict":
            pass

        else:
            print("'type' must be 'random' or 'dict'")


if __name__ == "__main__":
    # generate random chromosome
	init_chromosome =  Build()
	Build.random(init_chromosome)