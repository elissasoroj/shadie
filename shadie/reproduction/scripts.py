#!/usr/bin/env python

"""
Generic scripts.
"""
#---------------------------
#early ()
#standard early script
EARLY = """
// diploids (p1) just generated haploid gametophytes
    if (sim.generation % 2 == 0) {{

        // fitness affects gametophyte survival
        p0.fitnessScaling = (gam_pop_size / p0.individualCount);

        //set mutation rate for haploids
        sim.chromosome.setMutationRate(gam_mutation_rate);

        // p0 and p1 survival callbacks
        s1.active = 1;
        s2.active = 0;
        s3.active = 1;
        s4.active = 0;

        // haploids get modified fitness, without dominance
        {activate}
    }}


    // odd generations = gametophytes (p0) just generated sporophytes
    else {{

        // fitness affects sporophytes
        p1.fitnessScaling = spo_pop_size / p1.individualCount;

        //set mutation rate for diploids
        sim.chromosome.setMutationRate(spo_mutation_rate);

        // turn off p0 survival callbacks
        // turn on p1 survival callbacks
        s1.active = 0;
        s2.active = 1;
        s3.active = 0;
        s4.active = 1;

        // diploids get SLiM's standard fitness calculation, with dominance
        {deactivate}
    }}
"""

# gam_pop_size
# gam_female_to_male_ratio
# gam_mutation_rate
# gam_sperm_per_microscore
# sperm_pool (?? TODO)
# spo_pops_size
# spo_mutation_rate

EARLY_OPT = """
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

#-----------------------------------------------
#activate/deactivate fitness callbacks

ACTIVATE = "{idx}.active = 1;"

DEACTIVATE = "{idx}.active = 0;"

#-----------------------------------------------
#proper substitution behavior for gametophyte generation -
#must be repreated for every mut
SUBSTITUTION = """
    // gametophytes have just undergone fitness selection
    if (sim.generation % 2 == 0) {{
        {muts}
    }}
"""

SUB_MUTS = """
        mut{idx} = sim.mutationsOfType({mut});
        freq{idx} = sim.mutationFrequencies(NULL, mut{idx});
        if (any(freq{idx} == 0.5))
            sim.subpopulations.genomes.removeMutations(mut{idx}[freq{idx} == 0.5], T);
"""

#-----------------------------------------------
# gam_maternal_effect
GAM_MATERNAL_EFFECT_ON_P1 = """
	// Gametophyte mother fitness affects sporophyte survival
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * gam_maternal_effect) + fitness * (1 - gam_maternal_effect);
        return (draw < corrected_fitness);
    }
"""


# spo_maternal_effect
SPO_MATERNAL_EFFECT_ON_P0 = """
    // Sporophyte motherr fitness affects gametophyte survival
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * spo_maternal_effect) + fitness * (1 - spo_maternal_effect);
        return (draw < corrected_fitness);
    }
"""

#-----------------------------------------------
# spo_random_death_chance
# gam_random_death_chance

SURV = """
// remove p1 individuals during even generations
s1 survival(p1) {{
    return F;
}}

// remove p1s by random chance of death and apply maternal effects to fitness
s2 survival(p1) {{
    if (runif(1) < spo_random_death_chance)
        return F;
    {p1_maternal_effect}
    return NULL;
}}

// remove p0 individuals during odd generations
s3 survival(p0) {{
    return F;
}}

// remove p0s by random chance of death and apply maternal effects to fitness
s4 survival(p0) {{
    //this code implements random chance of death in gametophytes
    if (runif(1) < gam_random_death_chance)
        return F;   
    {p0survival} 
    {p0_maternal_effect}
    return NULL;
}}
"""