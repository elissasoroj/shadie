#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import 

class Phylogeny:
    def __init__(self, tree):
        
        # store input params
        self.tree = tree
        
        # will be used to store output results
        self.demo = #demography in proper format


# define a demographic model and the total genome size and recomb rate or map
mod = Model(
    tree=tree, 
    Ne=1000, 
    nsamples=2,          # number of sampled haplotypes per tip in final data 
    organism="fern",    # defines how gametes get selected and replicate
    recomb=1e-9, 
    genome_size=1e6,
)

# call simulate with details on genome structure (and which life stage selection occurs?)
mod.simulate(
    dict(name='a', selection=-0.01, start=2000, end=3000, haploid=True),
    dict(name='b', selection=-0.02, start=5000, end=6000, haploid=True),
    dict(name='c', selection=-0.03, start=9000, end=10000, haploid=True),
    dict(name='A', selection=0.01, start=2000, end=3000, haploid=False),
    dict(name='B', selection=0.02, start=5000, end=6000, haploid=False),
    dict(name='C', selection=0.03, start=9000, end=10000, haploid=False),
)

#demography code from ipcoal:

# demography info to fill
        self.ms_migrate = []
        self.ms_migtime = []
        self.ms_demography = set()
        self.ms_popconfig = ms.PopulationConfiguration()

        # get migration time, rate {mrates: [], mtimes: []}
        self._get_migration()

        # get .ms_demography dict for msprime input
        self._get_demography()

        # get .ms_popconfig as msprime input
        self._get_popconfig()

        # this is used when tips are not ultrametric
        self._get_nsamples()

        # to hold the model outputs
        self.df = None
        self.seqs = None
        self.ancestral_seq = None

        # check substitution model kwargs and assert it is a dict
        self.substitution_model = substitution_model
        self._check_substitution_kwargs()