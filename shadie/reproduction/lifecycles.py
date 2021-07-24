#!/usr/bin/env python

"""
Convenience functions for constructing each reproduction mode
"""

from shadie.reproduction.scripts import FIT, ACTIVATE, DEACTIVATE, EARLY
#from shadie.chromosome.build import ChromosomeRandom, Chromosome, ChromosomeExplicit


class Bryophyte:
    def __init__(self):
        return None

    def dioicous(
        self,
        chromosome,
        dNe = int(1000), #diploid Ne
        gNe =  int(2000),  #haploid Ne
        ftom = float.as_integer_ratio(1/1), #F:M sporophyte
        spores = int(100), #number of spores per sporopyte
        clonerate = float(1.0), #chance of cloning
        clones = float(1.0), #number of clones
        selfrate = float(0.0), #chance of intragametophytic selfing
        maternalweight = float(0.0), #maternal contribution to fitness
        deathchance = float(0.0), #random chance of death for both stages
        _maletag = 0,
        _femtag = 1, 
        _usedtag = 2,
        ):

        self.chromosome = chromosome

        #fitness callback:
        i = 4
        for mut in self.chromosome.mutations:
            i = i + 1
            idx = str("s"+str(i))
            self.fitdict = {(mut, idx), FIT}

        for mut in self.fitdict:
            self.activdict = {mut[1], ACTIVATE}
            self.deactivdict = {mut[1], DEACTIVATE}

        #activedict
        for i in self.activdict:
            activate_str = "\n  ".join([i.strip(';') + ';'])

        for i in self.deactivdict:
            deactivate_str = "\n  ".join([i.strip(';') + ';'])

        early_script = (
            EARLY.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        return([None, early_script]) #????? 


    def monoicous(
        self,
        chromosome,
        dNe = int(1000), #diploid Ne
        gNe =  int(2000),  #haploid Ne
        ftom = float.as_integer_ratio(1/1), #F:M sporophyte
        spores = int(100), #number of spores per sporopyte
        clonerate = float(1.0), #chance of cloning
        clones = float(1.0), #number of clones
        selfrate = float(0.0), #chance of intragametophytic selfing
        maternalweight = float(0.0), #maternal contribution to fitness
        deathchance = float(0.0), #random chance of death for both stages
        ):
        """
        ...
        """
        return("")


class Pteridophyte:
    def __init__(self):
        return None
    """
    Reproduction mode based on ferns and lycophytes
    """

    def moneicous(self):
        """
        ...
        """
        return("")

class Spermatophyte:
    def __init__(self):
        return None


if __name__ == "__main__":
    import shadie

    # define mutation types
    m0 = shadie.mtype(0.5, 'n', 2.0, 1.0)
    m1 = shadie.mtype(0.5, 'g', 3.0, 1.0)
    m2 = shadie.mtype(0.5, 'f', 0)
    
    # define elements types
    e0 = shadie.etype([m0, m1], [1, 2])
    e1 = shadie.etype([m2], [1])

    # design chromosome of elements
    # Do we want users to be able to put in a chromosome like this 
    #and have the gaps filled with neutral portions?
    chrom = shadie.chromosome.explicit({
        (500, 1000): e1,
        (2000, 3000): e0,
        (3001, 5000): e1,
    })

    repro = shadie.reproduction.Bryophyte().dioicous(chromosome = chrom)
    print(repro)