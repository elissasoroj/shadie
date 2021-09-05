#!/usr/bin/env python

"""
Generic scripts.
"""

#
# gam_pop_size
# gam_female_to_male_ratio
# gam_mutation_rate
# gam_sperm_per_microscore
# sperm_pool (?? TODO)
# spo_pops_size
# spo_mutation_rate
# 
EARLY = """
	// even generation, sporophytes (p1) just generated gametophytes
    if (sim.generation % 2 == 0) {{
        
        // set fitness to affect (female) gametophytes
        megaspores = p0.individuals[p0.individuals.tag==1];
        megaspores.fitnessScaling = gam_pop_size * gam_female_to_male_ratio / length(megaspores);
        
        // set fitness to affect male gametophytes
        // extra pressure applied to sperm to reduce sim size
        microspores = p0.individuals[p0.individuals.tag==2];
        microspore_pool = 2 * sperm_pool / gam_sperm_per_microspore;
        microspores.fitnessScaling = microspore_pool / length(microspores);

        // set mutation rate for haploids
        sim.chromosome.setMutationRate(gam_mutation_rate);

        // set survival callbacks off for p1 and on for p0
        s1.active = 1;
        s2.active = 0;
        s3.active = 1;
        s4.active = 0;

        // set modified fitness calc (w/o dominance) to gametophytes for each MutationType.
        {activate}
    }}

    // odd generation, gametophytes (p0) just generated sporophytes
    else {{

        // set fitness to affect sporophytes
        p1.fitnessScaling = spo_pop_size / p1.individualCount;

        // set mutation rate to sporophyte rate
        sim.chromosome.setMutationRate(spo_mutation_rate);
        
        // set survival callbacks off for p0 and on for p1
        s1.active = 0;
        s2.active = 1;
        s3.active = 0;
        s4.active = 1;

        // set standard fitness calc (w/ dominance) to sporophytes for each MutationType
        {deactivate}
    }}
"""

# Is this not redundant with SURV?
SPO_DEATH_CHANCE = """
// apply random chance of death in sporophytes
    if (runif(1) < spo_random_death_chance)
        return F;
    else
        return NULL;
"""

# Is this not redundant with SURV?
GAM_DEATH_CHANCE = """
// apply random chance of death in gametophytes
    if (runif(1) < gam_random_death_chance)
        return F;
    else
        return NULL;
"""


#
# gam_maternal_effect
#
MATERNAL_EFFECT_P0 = """
	// get maternal effect as weighted average
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * gam_maternal_effect) + fitness * (1 - gam_maternal_effect);
        return (draw < corrected_fitness);
    }
"""

#
# spo_maternal_effect
#
MATERNAL_EFFECT_P1 = """
    // maternal effect as weighted average
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * spo_maternal_effect) + fitness * (1 - spo_maternal_effect);
        return (draw < corrected_fitness);
    }
"""


#
# spo_random_death_chance
# gam_random_death_chance
# 
#
SURV = """
// remove p1 individuals during even generations
s1 survival(p1) {{
    return F;
}}

// remove p1s by random chance of death and apply maternal effects to fitness
s2 survival(p1) {{
    if (runif(1) < spo_random_death_chance)
        return F;
    else
        return NULL;

    {p0_maternal_effect}
}}

// remove p0s by random chance of death and apply maternal effects to fitness
s3 survival(p0) {{
    if (runif(1) < gam_random_death_chance)
        return F;
    else
        {p0_survival}
    
    {p1_maternal_effect}
}}

// remove p0 individuals during odd generations
s4 survival(p0) {{
    return F;
}}
"""

ANGIO_SURV_P0 = """
    {//All ovules survive; this is a way of implementing maternal effect
    //if mother died, they would not be produced
        if (individual.tag == 1)
                return T;
        else
            return NULL;
    }
"""