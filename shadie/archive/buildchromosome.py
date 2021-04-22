#!/usr/bin/env python

"""
Generates chromosome strurcture for SLiM simulation
"""

from typing import Union
import pandas as pd
import numpy as np
from loguru import logger

# internal imports
from shadie.mutations import MutationList
from shadie.elements import ElementList
from shadie.elements import ElementType

# default element types
from shadie.globals import NONCDS, INTRON, EXON
from shadie.globals import NEUT, SYN, DEL, BEN


class ChromosomeBase:
    def __init__(self, genome_size):
        self.genome_size = genome_size
        self.data = pd.DataFrame(
            columns=['name', 'start', 'end', 'eltype', 'script'],
            data=None,
        )       


class ChromosomeRandom(ChromosomeBase): 
    """
    Generates a random chromosome from of a given length from a set of
    intron, exon, and non-cds genomic ElementType objects. The default
    elements are used if not entered by the user.
    """
    def __init__(
        self, 
        genome_size:int=20000, 
        intron:ElementType=None,
        exon:ElementType=None,
        noncds:ElementType=None,
        seed:Union[int, None]=None,
        ):

        super().__init__(genome_size)
        self.rng = np.random.default_rng(seed)
        self.intron = intron if intron is not None else INTRON
        self.exon = exon if exon is not None else EXON
        self.noncds = noncds if noncds is not None else NONCDS
        self.run()


    def get_noncds_span(self, scale:int=5000):
        """
        Draws the number of bases until the next element from an 
        exponential distribution. The scale is the average waiting
        time in number of bp.
        """
        return int(self.rng.exponential(scale=scale))


    def get_cds_spans(self, length_scale:int=1000, intron_scale:int=1000):
        """
        Draws the number of exons in a fixed length space from a 
        Poisson distribution. The lam parameter is the average number
        of events per sampled region. A value of 0.005 means one intron
        per 200bp.
        """
        cds_span = int(self.rng.exponential(scale=length_scale))
        n_introns = int(self.rng.poisson(cds_span / intron_scale))
        if n_introns:
            splits = self.rng.dirichlet(np.ones(n_introns * 2 - 1))
            splits = (splits * cds_span).astype(int)
            splits[-1] = cds_span - sum(splits[:-1])
        else:
            splits = np.array([cds_span])
        return splits


    def run(self, noncds_scale=5000, cds_scale=1000, intron_scale=1000):
        """
        Generates a chromosome by randomly sampling waiting times 
        between CDS regions, and the number of introns within CDS
        regions.
        """
        idx = 0
        while 1:
            # get non-cds span
            pos = self.get_noncds_span(noncds_scale)
            self.data.loc[idx] = (
                self.noncds.altname, 
                idx, min(idx + pos, self.genome_size), 
                self.noncds.name, self.noncds,
            )
            idx += pos
            
            # get cds span
            posses = self.get_cds_spans(cds_scale, intron_scale)

            # break if cds goes beyond the end of the genome.
            if idx + posses.sum() > self.genome_size:
                break

            # enter the cds into data
            for enum, pos in enumerate(posses):
                if not enum % 2:
                    self.data.loc[idx] = (
                        self.exon.altname, 
                        idx, idx + pos, 
                        self.exon.name, self.exon,
                    )
                else:
                    self.data.loc[idx] = (
                        self.intron.altname,
                        idx, idx + pos, 
                        self.intron.name, self.intron,
                    )
                idx += pos



class ChromosomeStandard(ChromosomeBase):
    """
    Builds the default shadie chromosome used for testing.
    """
    def __init__(self):
        super().__init__(genome_size=100000)
        self.elements = ElementList(NONCDS, INTRON, EXON)
        self.run()

    def run(self):
        """
        Generates a chromosome of a fixed length and composition.
        """
        idx = 0
        while 1:
            self.data.iloc[idx] = (
                NONCDS.altname, idx, idx + 2000, NONCDS.name, NONCDS)
            idx += 2000
            self.data.iloc[idx] = (
                EXON.altname, idx, idx + 2000, EXON.name, EXON)
            idx += 2000
            self.data.iloc[idx] = (
                INTRON.altname, idx, idx + 2000, INTRON.name, INTRON)
            idx += 2000
            self.data.iloc[idx] = (
                EXON.altname, idx, idx + 2000, EXON.name, EXON)
            idx += 2000
            self.data.iloc[idx] = (
                NONCDS.altname, idx, idx + 2000, NONCDS.name, NONCDS)
            if idx >= self.genome_size:
                break


class ChromosomeExplicit(ChromosomeBase):
    """
    Builds a chromosome dataframe from explicit instructions provided
    as start, stop positions of ElementTypes.

    chromosome.explicit({
        500: g1,
        1000: g2,
        1200: g1,
        1500: g3,
        2000: g2,
    })
    """
    def __init__(self, data):
        super().__init__(genome_size=max(data.keys()))

        # check data dict for proper structure
        assert all(isinstance(i, int) for i in data.keys()), (
            "keys of input data should be integers.")
        assert all(isinstance(i, ElementType) for i in data.values()), (
            "values of input data should be ElementType objects.")

        # entere explicit dict into data
        idx = 0
        for key in sorted(data):
            self.data.loc[idx] = (
                data[key].altname, idx, idx + key, data[key].name, data[key])



class Build:
    """
    Chromosome object created based on user parameters

    Parameters:
    -----------
    chromtype(str): default = "random"
        "random" will generate a random chromosome
        "dict" will accept a dictionary object

    exoncount(int): default = None
        Will produce a random chromosome using defaults with specified number 
        of genes (composed of exons only)

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

    exons (ElementList): default = None
       list of genomic elements that may be used as exons

    introns (ElementList): default = None
       list of genomic elements that may be used as introns

    noncoding (ElementList): default = None
       list of genomic elements that may be used as non-coding regions

    mutation_rate(int): default = 1e-7
        Chance of mutation at each nucleotide

    genome_size(int): default = =1e6
        Length of chromosome in nucleootides
    """    
    def __init__(
    	self,
    	chromtype="random",  #"random" or "dict" 
        genecount=None,
		exons=None,          # list of ElementTypes eligible for exon regions
        introns=None,        # list of ElementTypes eligible for intron regions
        noncoding=None,      # list of ElementTypes eligible for non-coding regions
		elementlist=None,    # read in ElementList object
		genome_size=1e6,     # will be used to calculate chromosome end (length -1)     
        ):
    
        self.type = chromtype
        self.exoncount = genecount
        self.exons = exons
        self.introns = introns
        self.noncoding = noncoding
        self.elements = elementlist
        self.gensize = genome_size

        self.genedf = pd.DataFrame([])



        # random generation attrs
        self._maxexon = int(self.gensize / (self.exoncount * 1.25))
        self._minexon = int(self.gensize / (self.exoncount * 3))
        self._ncmin = int((self.gensize - (self.exoncount * maxexon)) / (self.exoncount + 1))
        self._ncmax = int((self.gensize - (self.exoncount * minexon)) / (self.exoncount + 1))


        # if self.type == "custom":
        #     if isinstance(args, io.TextIOBase):
        #         self.genedf = pd.read_csv(args)

        #     elif isinstance(args, pd.DataFrame):
        #         self.genedf = args

        #     else:
        #         genelist = []
        #         for i in args:
        #             if isinstance(i, dict):
        #                 genelist.append(i)
        #         self.genedf = pd.DataFrame(genelist)

        # elif self.type == "random":
        #     pass
        

    def dict(self, genedf):
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


    def generate_genes(self):
        """
        Generates a chromosome with specified number of exons
        """
        # exon size

        while 1:
            exons = 0
            idx = 0

            # make initial non-coding region
            nc_length = np.random.randint(ncmin, ncmax)
            self.df.loc = (
                "noncoding",
                None,
                idx,
                idx + nc_length + 1,
                NONCOD.name,
                NONCOD,
            )
            idx = idx + nc_length
                
            # make first exon
            ex_length = np.random.randint(minexon, maxexon) + 1
            self.df.loc[idx] = (
                "exon", 
                None, 
                idx, 
                idx + ex_length - 1, 
                EXON.name, 
                EXON,
            )
            idx = idx + ex_length
            exons += 1

            # add additional exons to reach sampled count
            while exons < self.exoncount:
                nc_length = np.random.randint(ncmin, ncmax)
                genelements.loc[base, 'type'] = "noncoding"
                genelements.loc[base, 'eltype'] = NONCOD.name
                genelements.loc[base, 'script'] = NONCOD
                genelements.loc[base, 'start'] = base
                genelements.loc[base, 'finish'] = base + nc_length - 1
                base = base + nc_length

                #make first exon
                ex_length = np.random.randint(minexon, maxexon) + 1
                genelements.loc[base, 'type'] = "exon"
                genelements.loc[base, 'eltype'] = EXON.name
                genelements.loc[base, 'script'] = EXON
                genelements.loc[base, 'start'] = base
                genelements.loc[base, 'finish'] = base + ex_length -1
                base = base + ex_length
                exons += 1

                #final non-coding region
                genelements.loc[base, 'type'] = "noncoding"
                genelements.loc[base, 'eltype'] = NONCOD.name
                genelements.loc[base, 'script'] = NONCOD
                genelements.loc[base, 'start'] = base
                genelements.loc[base, 'finish'] = self.gensize - 1

                if base > (self.gensize-1):
                    logger.info("trying again")
                    exons = 0
                    base = int(0)
                    continue
                else:
                    break

                logger.debug(genelements)
                self.genelements = genelements
                self.mutationlist = MutationList(NEUT, SYN, DEL, BEN)
                self.elementlist = ElementList(self.mutationlist, EXON, NONCOD)



    def random(self):
        "generates a random chromosome"

        if self.type == "random":
            if self.exons is None:
                genelements = pd.DataFrame(
                columns=['type', 'name', 'start', 'finish', 'eltype', 'script'],
                data=None,
                )
                finalnc_length = np.random.randint(3000, 5000)
                end = self.gensize - finalnc_length
                logger.debug("Made objects: base = {base}, genelements = {}, end = {end}")

                while True:
                    base = int(0)
                    while base < end:
                        #make initial non-coding region
                        nc_length = np.random.randint(100, 5000)
                        genelements.loc[base, 'type'] = "noncoding"
                        genelements.loc[base, 'eltype'] = NONCOD.name
                        genelements.loc[base, 'script'] = NONCOD
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + nc_length - 1
                        base = base + nc_length
                    
                        #make first exon
                        ex_length = round(np.random.lognormal(np.log(250), np.log(1.3))) + 1
                        genelements.loc[base, 'type'] = "exon"
                        genelements.loc[base, 'eltype'] = EXON.name
                        genelements.loc[base, 'script'] = EXON
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + ex_length -1
                        base = base + ex_length
                        logger.info("Gene added")    
                        
                        while np.random.random_sample() < 0.75:  #25% probability of stopping
                            in_length = round(np.random.normal(450, 100))
                            genelements.loc[base, 'type'] = "intron"
                            genelements.loc[base, 'eltype'] = INTRON.name
                            genelements.loc[base, 'script'] = INTRON
                            genelements.loc[base, 'start'] = base
                            genelements.loc[base, 'finish'] = base + in_length -1
                            base = base + in_length
                          
                            ex_length = round(np.random.lognormal(np.log(250), np.log(1.3))) + 1
                            genelements.loc[base, 'type'] = "exon"
                            genelements.loc[base, 'eltype'] = EXON.name
                            genelements.loc[base, 'script'] = EXON
                            genelements.loc[base, 'start'] = base
                            genelements.loc[base, 'finish'] = base + ex_length -1
                            base = base + ex_length 
                              
                    #final non-coding region
                    genelements.loc[base, 'type'] = "noncoding"
                    genelements.loc[base, 'eltype'] = NONCOD.name
                    genelements.loc[base, 'script'] = NONCOD
                    genelements.loc[base, 'start'] = base
                    genelements.loc[base, 'finish'] = self.gensize - 1

                    logger.debug(genelements)
                    self.genelements = genelements
                    self.mutationlist = MutationList(NEUT, SYN, DEL, BEN)
                    self.elementlist = ElementList(self.mutationlist, EXON, INTRON, NONCOD)

                    if base > (self.gensize-1):
                        logger.info("trying again")
                        base = int(0)
                        continue
                    else:
                        break
       
            elif self.exons is not None: 
                #check the types
                for i in self.exons:
                    if isinstance(i, ElementType):
                        pass
                    else:
                        raise TypeError("exons must be ElementType class objects")

                for i in self.introns:
                    if isinstance(i, ElementType):
                        pass
                    else:
                        raise TypeError("introns must be ElementType class objects")

                #check the types
                for i in self.noncoding:
                    if isinstance(i, ElementType):
                        pass
                    else:
                        raise TypeError("noncoding must be ElementType class objects")

                genelements = pd.DataFrame(
                columns=['name', 'start', 'finish', 'eltype', 'script'],
                data=None,
                )
                finalnc_length = np.random.randint(self.gensize/300, self.gensize/200)
                end = self.gensize - finalnc_length
                logger.debug("Made objects: base = {base}, genelements = {}, end = {end}")

                while True:
                    base = int(0)
                    while base < end:
                        #make initial non-coding region
                        nc_length = np.random.randint(100, 5000)
                        choose = random.choice(self.noncoding)

                        genelements.loc[base, 'type'] = "noncoding"
                        genelements.loc[base, 'name'] = choose.altname
                        genelements.loc[base, 'eltype'] = choose.name
                        genelements.loc[base, 'script'] = self.elements.elementdict[choose.name]
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + nc_length - 1
                        base = base + nc_length

                        #make first exon
                        ex_length = round(np.random.lognormal(np.log(250), np.log(1.3))) + 1
                        choose = random.choice(self.exons)

                        genelements.loc[base, 'type'] = "exon"
                        genelements.loc[base, 'name'] = choose.altname
                        genelements.loc[base, 'eltype'] = choose.name
                        genelements.loc[base, 'script'] = self.elements.elementdict[choose.name]
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + ex_length -1
                        base = base + ex_length
                        logger.info("Gene added")    

                        while np.random.random_sample() < 0.75:  #25% probability of stopping
                            in_length = round(np.random.normal(450, 100))
                            choose = random.choice(self.introns)
                            genelements.loc[base, 'type'] = "intron"
                            genelements.loc[base, 'name'] = choose.altname
                            genelements.loc[base, 'eltype'] = choose.name
                            genelements.loc[base, 'script'] = self.elements.elementdict[choose.name]
                            genelements.loc[base, 'start'] = base
                            genelements.loc[base, 'finish'] = base + in_length -1
                            base = base + in_length

                            ex_length = round(np.random.lognormal(np.log(250), np.log(1.3))) + 1
                            choose = random.choice(self.exons)

                            genelements.loc[base, 'type'] = "exon"
                            genelements.loc[base, 'name'] = choose.altname
                            genelements.loc[base, 'eltype'] = choose.name
                            genelements.loc[base, 'script'] = self.elements.elementdict[choose.name]
                            genelements.loc[base, 'start'] = base
                            genelements.loc[base, 'finish'] = base + ex_length -1
                            base = base + ex_length 

                    #final non-coding region
                    choose = random.choice(self.noncoding)

                    genelements.loc[base, 'type'] = "noncoding"
                    genelements.loc[base, 'name'] = choose.altname
                    genelements.loc[base, 'eltype'] = choose.name
                    genelements.loc[base, 'script'] = self.elements.elementdict[choose.name]
                    genelements.loc[base, 'start'] = base
                    genelements.loc[base, 'finish'] = self.gensize - 1
                    
                    if base > (self.gensize-1):
                        logger.info("trying again")
                        base = int(0)
                        continue
                    else:
                        break

                logger.debug(genelements)
                self.genelements = genelements
                self.mutationlist = self.elements.mutationlist
                self.elementlist = self.elements

        else:
            print("'type' must be 'random' or 'dict'")


if __name__ == "__main__":

    #test custom builder:
    from shadie.mutations import MutationType

    mut1 = MutationType(0.5, "f", .03)
    mut2 = MutationType(0.5, "e", 0.4)
    mutlist = MutationList(mut1, mut2)
    eltype1 = ElementType(mut1, 1)
    eltype2 = ElementType(mut2, 1)
    mylist = ElementList(mutlist, eltype1, eltype2)

    custom = Build(
        exons = [eltype1, eltype2], 
        introns = [eltype1, eltype2], 
        noncoding = [eltype1, eltype2],
        elementlist = mylist,
    )

    # generate random chromosome
    init_chromosome =  Build()
    Build.random(custom)

    #chromosome with specified number of genes
    fourgenes = Build(genecount = 4)
    fourgenes.genes()
