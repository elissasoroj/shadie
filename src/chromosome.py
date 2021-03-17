#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import pandas as pd

class Chromosome:
	"""
    Chromosome object created based on user parameters
	"""
	#standard mutation types, with structure:
	#(name, dominance, distribution, {following depends on distribution})
	NEUT = ("m100", 0.5, "f", 0.0)		#neutral mutation
	SYN = ("m101", 0.5, "f", 0.0)			#synonymous
	DEL = ("m102", 0.1, "g", -0.03, 0.2)	#deleterious
	BEN = ("m103", 0.8, "e", 0.1)			#beneficial

	#standard genetic elements have the structure:
	#(name, mutations, frequency of each mutation)
	EXON = ("g100", c(m100,m101,m102), c(2,8,0.1))	#exon
	INTRON = ("g101", c(m100,m102), c(9,1))			#intron
	NONCOD = ("g102", c(m100), 1) 				#non-coding

	def __init__(
		self,
		genes=1,
		introns = 0,
		exons = 1,
		):
	"""
	Builds the chromosome for SLiM3 simulation

	Parameters:
	-----------
	genes(int): default = 1
		Number of genes on chromosome
		Can be supplied as `int` for random chromosome generation or
		As `dict` for explicitly defines genes, e.g:
		# call simulate with details on genome structure (and which life stage selection occurs?)
mod.simulate(
    dict(name='a', selection=-0.01, start=2000, end=3000, haploid=True),
    dict(name='b', selection=-0.02, start=5000, end=6000, haploid=True),
    dict(name='c', selection=-0.03, start=9000, end=10000, haploid=True),
    dict(name='A', selection=0.01, start=2000, end=3000, haploid=False),
    dict(name='B', selection=0.02, start=5000, end=6000, haploid=False),
    dict(name='C', selection=0.03, start=9000, end=10000, haploid=False),
)

	dict(name = "name", mutations = (fitness effect,  dominance), start = bp, end = bp)
	dict(name = "a", mutations = (-.01, 0.5), start = 2000, end = 3000)

	introns (int): default = 0
		Number of introns per gene

	exons (int): default = 1
		Number of exons
	"""

	def make():
	if isinstance(self.genes, int):
		print("genes must be dict object")

	elif isinstance(self.genes, dict):
		self.genedict = self.genes
		for self.gene in self.genedict:
			#make the mutation types
			mutdict = {}
			#...
			self.mutdict = mutdict

			#make the genomic element types
			elemdict = {}
			#...
			self.selemdict = elemdict

			#write the initializeGenomicElement lines to a dict
			initdict = {}
			#...
			self.initdict = initdict


	def make_rand():
	if isinstance(self.genes, int):
		pass

	def review():
		"prints mutations, genetic elements, and genes dataframe"
		#mutations
		print("Mutations:\n" self.mutdict)
		#genelems
		#genes (pandas dataframe)

if __name__ == "__main__":

    # make a random tree with 10 tips and root height 1M
    chrom1 = Chromosome(genes = 5)
    Chromosome.make(chrom1)
  
