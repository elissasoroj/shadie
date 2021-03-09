#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import pandas as pd


class Shadie(object):
	"""
    A program.script object for writing the script for simulation in SLiM3
	"""
	def __init__(
		self,
		tree=None,
		Ne = 1000,			#K
		nsamples=2,         # number of sampled haplotypes per tip in final data 
		organism="pter",    # defines how gametes get selected and replicate
		mutrate=1e-7,
		recomb=1e-9,		#sets rate in `initializeRecombinationRate`, also accepts map
		genome_size=1e6,	#will be used to calculate chromosome end (length -1)
		model = "nonWF",	#write "nonWF" to initailizeSLiMModelType()
		treeseq = "T", #turns on tree sequence recording 
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

		....
	"""

	def simulate(self):
		"calls SLiM to run the simulations"

	def write(self, filename="shadie.slim"):
		"writes the .slim script; optional to provide filename as 'filename.slim'"
	    filename.self = filename
	    mutrate.self = mutrate
	    muttype.self = muttype
	    geneltype.self = geneltype
	    genel.self = genel
	    recomb.self = recomb
	    Ne.self = Ne
	    model.self = model
	    #write initialize callbacks
	    script = open(filename.self, "a") #appends so that user does not accidentally overwrite old simulation
	    L1 = (
	    	"initialize() {\ninitializeSLiMModelType("+model.self+");\n"
	    	f"defineConstant('K',{Ne.self});\n"
	    	f"initializeMutationRate({mutrate.self});\n"
	    	f"initializeMutationType{muttype.self};\n m1.convertToSubstitution = T;\n"
	    	f"initializeGenomicElementType{geneltype.self};\n"
	    	f"initializeGenomicElement{genel.self};\n"
	    	"initializeRecombinationRate("+recomb.self+");\n}"
	    )
	    script.write(L1)
	    script.close

	def genome(self):
		if genes = 1:



# call simulate with details on genome structure (and which life stage selection occurs?)
mod.simulate(
    dict(name='a', selection=-0.01, start=2000, end=3000, haploid=True),
    dict(name='b', selection=-0.02, start=5000, end=6000, haploid=True),
    dict(name='c', selection=-0.03, start=9000, end=10000, haploid=True),
    dict(name='A', selection=0.01, start=2000, end=3000, haploid=False),
    dict(name='B', selection=0.02, start=5000, end=6000, haploid=False),
    dict(name='C', selection=0.03, start=9000, end=10000, haploid=False),
)


	def organism(self):
		"defines haploid/diploid life cycle offspring generation, selection, fitness effects"

		if self.organism = "pter":
			#define reproduction - sporophyte produces male or hermaphroditic gametophytes
			#define sex ratio?  - this may be relevant (e.g. hermaphroodite gametophytes 
				#induce formation of male-only gametophyte around them)
			#define fitness effects - selection on gametophyte
			#**no selfing with male/female??

		elif self.organism = "bryo":
			#define reproduction - male (sporophyte only?), female
			#define sex ratio?  - females live longer than males, 
			#define fitness effects - selection on gametophyte occurrs for longer than selection 
				# on sporophyte
			#**no selfing with male/female??

		elif self.organism = "angio":
			#define reproduction - hermaphrodite spermatophytes, male/female gameotphyte
			#define sex ratio?  - can sex ratio be set for haploid? male >> female for gametophyte
			#define fitness effects - selection mainly acts on sporophyte
			

		else:
			#base SLiM hermaphrodite? 

class Genome(object):
	"""
    Genome object created based on user parameters
	"""
	#standard mutation types
	NEUT = ("m1", 0.5, "f", 0.0)	#neutral mutation
	SYN = ("m2", 0.5, "f", 0.0)		#synonymous
	DEL = ("m3", 0.1, "g", -0.03, 0.2)	#deleterious
	BEN = ("m4", 0.8, "e", 0.1)		#beneficial

	#standard genetic elements
	EXON = ("g1", c(m2,m3,m4), c(2,8,0.1))	#exon
	INTRON = ("g2", c(m1,m3), c(9,1))		#intron
	NONC = ("g3", c(m1), 1) 				#non-coding

	def __init__(
		self,
		genes=1,
		introns = 0,
		exons = 0,

		):
	"""
	Builds the chromosome for SLiM3 simulation

	Parameters:
	-----------
	genes(int): default = 1
		Number of genes on chromosome

	introns (int): default = 0
		Number of introns per gene

	exons (int): default = 0
		Number of exons
	"""

	def simulate(self):
		"calls SLiM to run the simulations"