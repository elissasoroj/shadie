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
        {p0_fitnessScaling};

        //set mutation rate for haploids
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);

        // p0 and p1 survival callbacks
        s1.active = 1;
        s2.active = 0;
        s3.active = 1;
        s4.active = 0;

        // haploids reproduce, diploids don't
        s5.active = 0;
        s6.active = 1;

        // haploids get modified fitness, without dominance
        {activate}
    }}


    // odd generations = gametophytes (p0) just generated sporophytes
    else {{

        // fitness affects sporophytes
        p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount;

        //set mutation rate for diploids
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);

        // turn off p0 survival callbacks
        // turn on p1 survival callbacks
        s1.active = 0;
        s2.active = 1;
        s3.active = 0;
        s4.active = 1;

        // diploids reproduce, haploids don't
        s5.active = 1;
        s6.active = 0;

        // diploids get SLiM's standard fitness calculation, with dominance
        {deactivate}
    }}
"""


P0_FITNESS_SCALE_DEFAULT = "p0.fitnessScaling = GAM_POP_SIZE / p0.individualCount;"

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
        megaspores.fitnessScaling = GAM_POP_SIZE * GAM_FEMALE_TO_MALE_RATIO / length(megaspores);
        
        // set fitness to affect male gametophytes
        // extra pressure applied to sperm to reduce sim size
        microspores = p0.individuals[p0.individuals.tag==2];
        microspore_pool = 2 * sperm_pool / GAM_SPERM_PER_MICROSPORE;
        microspores.fitnessScaling = microspore_pool / length(microspores);

        // set mutation rate for haploids
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);

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
        p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount;

        // set mutation rate to sporophyte rate
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);
        
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
    {late}
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
        corrected_fitness = (maternal_effect * GAM_MATERNAL_EFFECT) + fitness * (1 - GAM_MATERNAL_EFFECT);
        return (draw < corrected_fitness);
    }
"""


# spo_maternal_effect
SPO_MATERNAL_EFFECT_ON_P0 = """
    // Sporophyte motherr fitness affects gametophyte survival
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * SPO_MATERNAL_EFFECT) + fitness * (1 - SPO_MATERNAL_EFFECT);
        return (draw < corrected_fitness);
    }
"""

#-----------------------------------------------
# spo_random_death_chance
# gam_random_death_chance

SURV = """
// remove p1 individuals during even generations
s1 survival(p1) {{
    if ((individual.tag == 44) | (individual.tag == 5) | (individual.tag == 45)){{
        individual.tag = 0;
        return T;
    }}
    else
        return F;
}}

// remove p1s by random chance of death and apply maternal effects to fitness
s2 survival(p1) {{
    if (runif(1) < SPO_RANDOM_DEATH_CHANCE)
        return F;
    {p1_maternal_effect}
    return NULL;
}}

// remove p0s by random chance of death and apply maternal effects to fitness
s3 survival(p0) {{
    //this code implements random chance of death in gametophytes
    if (runif(1) < GAM_RANDOM_DEATH_CHANCE)
        return F;   
    {p0survival} 
    {p0_maternal_effect}
    return NULL;
}}

// remove p0 individuals during odd generations
s4 survival(p0) {{
    if ((individual.tag == 4) | (individual.tag == 6)) {{
        individual.tag = 0;
        return T;
    }}
    else
        return F;
}}
"""