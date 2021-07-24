#!/usr/bin/env python

"""
Convenience functions for constructing each reproduction mode
"""

from shadie.reproduction.scripts import FIT, ACTIVATE, DEACTIVATE, EARLY
#from shadie.chromosome.build import ChromosomeRandom, Chromosome, ChromosomeExplicit

REPRO_BRYO_DIO_p1 = """
{{  // creation of spores from sporophytes
    g_1 = genome1;
    g_2 = genome2;
    
    meiosis_reps = floor(Spore_num/2);
    for (rep in 1:meiosis_reps)
    {{
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = ifelse (runif(1)<FtoM, 1, 0);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = ifelse (runif(1)<FtoM, 1, 0);
    }}

}}
"""

REPRO_BRYO_DIO_p0 = """
{{  //creation of sporophyte from haploids
    if (individual.tag == 1)    // females find male gametes to reproduce
    {{
        reproduction_opportunity_count = 1;
        
        // clones give the focal individual extra opportunities to reproduce
        if (runif(1) <= Clone_rate)
            reproduction_opportunity_count = reproduction_opportunity_count + Clone_num;
        
        for (repro in seqLen(reproduction_opportunity_count))
        {{
            if (runif(1) <= Self_rate)
            {{
                // this is selfing using two identical gametes – intragametophytic selfing
                // intergametophytic selfing might happen below, by chance
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
            }}
            else
            {{
                sperm = p0.sampleIndividuals(1, tag=0); // find a male!
                
                if (sperm.size() == 1)
                {{
                    child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                    
                    if (Maternal_weight > 0) //Mother's fitness affects sporophyte fitness; see survival()
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    
                    sperm.tag = 2;  // take out of the mating pool
                }}
            }}
        }}
    }}
}}
"""

REPRO_BRYO_MONO = """

"""

REPRO_ANGIO_DIO_p1 = """
{{ // creation of gametes from sporophytes
    g_1 = genome1;
    g_2 = genome2;
    
    if (individual.tag == 1)
    {
        // determine how many ovules were fertilized, out of the total
        fertilizedOvules = rbinom(1, ovule_count, fertilization_rate);
        meiosis_reps = floor(fertilizedOvules/2);
        if (runif(1) <= Clone_rate)
            meiosis_reps = meiosis_reps*Clone_num;
        
        for (rep in 1:meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
        }
    
    }
    else //individual is male
    {
        meiosis_reps = floor(pollen_count/2);
        if (runif(1) <= Clone_rate)
            meiosis_reps = meiosis_reps*2;
        for (rep in 1:meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
        }
    }
}}
"""

REPRO_ANGIO_DIO_p0  = """
{{ // creation of sporophytes from gametes
    if (individual.tag == 1)  // females find male gametes to reproduce
    {
        if (pollen_comp == T)
        {
            pollen_pool = p0.sampleIndividuals(pollen_per_stigma, tag=0);   // sperm land on stigma
            for (pollen in pollen_pool)
            {
                pollen.setValue("fitness", p0.cachedFitness(pollen.index)); //store fitness value
                pollen.tag = 2;
            }
            
            if (pollen_pool.length()>0)
            {
            target_fitness = max(pollen_pool.getValue("fitness"));
            winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
            sperm = winners[0];
            }
            else sperm = p0.sampleIndividuals(1, tag=0);    // find a male
        }
        else
            sperm = p0.sampleIndividuals(1, tag=0); // find a male
        if (sperm.size() == 1)
        {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
            sperm.tag = 2;
            
            if (runif(1) <= FtoM)
                child.tag = 1;
            else
                child.tag = 0;
        }
    }
}}
"""

REPRO_HOMOSPORE = """
//Hermaphroditic gametophytes
"""

class Bryophyte:
    def __init__(self):
        self.chromosome = None
        self.fitdict = None
        self.activdict = None
        self.deactivdict = None
        self.lifedict = None

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

        self.lifedict = {("early", None), early_script}
        self.lifedict = {("reproduction", "p1"), REPRO_BRYO_DIO_p1}
        self.lifedict = {("reproduction", "p0"), REPRO_BRYO_DIO_p0}


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

    repro = shadie.reproduction.Bryophyte()
    repro.dioicous(chromosome = chrom)
    print(repro.fitdict)