#!/usr/bin/env python

"""
Creates amino acid based script for SLiM

"""

#package imports


class Protein:
    """
    shadie.Chromosome object stores the genome information for a SLiM
    simulation. Allowds the use to inspect the genome parameters 
    before proceeding

    """

    def __init__(
        self,
        genome_size=2e3,        #will be used to calculate chromosome end (length -1)
        genome = None           #optional BuildChromosome object
        ):
        """
        Accepts a subclass Build object to define chromosome structure.
        Otherwise, takes genome_size (length) and mutation rate and 
        creates a single gene for the simulation using default EXON 
        genomic element type and default mutation types assigned to 
        EXON (see globals.py). 

        Parameters:
        -----------

        mutation_rate(int): default = 1e-7
            Chance of mutation at each bp


        genome_size(int): default = =1e6
            Length of chromosome
        """

        self.genome = genome

    def Translate(

    	)

    