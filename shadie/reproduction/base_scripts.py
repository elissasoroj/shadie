#!/usr/bin/env python

"""
String scripts for reproduction blocks
"""

EARLY_BRYO_DIO = """
    // even gens = diploids (p1) just generated haploid gametophytes
    if (sim.generation % 2 == 0) {
        
        // fitness affects gametophyte survival
        p0.fitnessScaling = (hK / p0.individualCount);
        
        // haploids get modified fitness, without dominance
        s1.active = 1;
        s2.active = 1;
        s3.active = 1;
        s4.active = 0;
        s5.active = 1;
        s6.active = 0;
    }

	// odd generations = gametophytes (p0) just generated sporophytes
    else {
        
        // fitness affects sporophytes
        p1.fitnessScaling = dK / p1.individualCount;

        // diploids get SLiM's standard fitness calculation, with dominance
        s1.active = 0;
        s2.active = 0;
        s3.active = 0;
        s4.active = 1;
        s5.active = 0;
        s6.active = 1;
    }
"""

FITNESS_BRYO_DIO_P1 = """
    // this code implements random death chance
    if (runif(1) < Death_chance)
        return F;
    else
        return NULL;
    
    // maternal effect
    maternal_effect = individual.getValue("maternal_fitness");
    
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * Maternal_weight) + fitness * (1 - Maternal_weight);
        return (draw < corrected_fitness);
    }
    return NULL;
"""  

FITNESS_BRYO_DIO_P0 = """
    // this code implements random death chance
    if (runif(1) < Death_chance)
        return F;
    else
        return NULL;
"""