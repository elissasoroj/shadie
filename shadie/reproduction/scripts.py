#!/usr/bin/env python

"""
Generic scripts.
"""

#standard first script

# ==PARAMETERS==
# SPO_MUTATION_RATE
# GAM_MUTATION_RATE
#-----------------
FIRST = """
// alternate generations
    if (sim.cycle % 2 == 1) {
        // reproduction(p0) will create SPOROPHYTES in p1.
        
        // set reproduction function to be used this generation     
        s0.active = 1;
        s1.active = 0;
        
        // set mutation rate that will apply to the offspring
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);
    }
    
    else {  //(sim.cycle % 2 == 0)
        // reproduction(p1) will produce GAMETOPHYTES into p0.
        
        // set reproduction function to be used this generation  
        s0.active = 0; //GAM
        s1.active = 1; //SPO
        
        // set mutation rate that will apply to the offspring
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);
    }
"""
#--------------------------------------------------------------
#standard early script

# ==PARAMETERS==
# GAM_RANDOM_DEATH_CHANCE
# GAM_CEILING
# GAM_POP_SIZE
# SPO_MATERNAL_EFFECT
# -------------------------
# TAGS
# 3
#-----------------
EARLY = """
// diploids (p1) just generated haploid gametophytes into p0
    if (community.tick % 2 == 0) {
        // alternate generations.
        
        //remove parents, except for clones
        non_clones = p1.individuals[p1.individuals.tag != 4];
        sim.killIndividuals(non_clones);
        //reset clone tags to regular tag
        p1.individuals.tag = 3;
        
        // new p0s removed by three possible mechanisms: 
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)
        
        // 1. make a mask for random death
        if (GAM_RANDOM_DEATH_CHANCE > 0) {
            random_death = (runif(p0.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p0.individuals[random_death]);
        }

        //random death also occurs to implement GAM_CEILING
        if (p0.individualCount > GAM_CEILING) {
            to_kill = GAM_CEILING - GAM_POP_SIZE;
            death_chance = to_kill/p0.individualCount;
            random_death = sample(c(F,T), p0.individualCount, T, c(1-death_chance, death_chance));
            sim.killIndividuals(p0.individuals[random_death]);
            }
        
        // 2. maternal effect re-calculation 
        if (SPO_MATERNAL_EFFECT > 0) {
            //temp fitness scaling
            scale = SPO_POP_SIZE / p1.individualCount;
            //calculations
            maternal_fitnesses = SPO_MATERNAL_EFFECT*p0.individuals.getValue("maternal_fitness");
            child_fitnesses = (1 - SPO_MATERNAL_EFFECT)*p0.cachedFitness(NULL);
            corrected_fitnesses = scale * (maternal_fitnesses + child_fitnesses);
            
            //remove inds based on maternal effects
            to_kill = (runif(p0.individualCount) > corrected_fitnesses);
            sim.killIndividuals(p0.individuals[to_kill]);
        }
        
        // fitness affects gametophyte survival
        p0.fitnessScaling = GAM_POP_SIZE / (p0.individualCount);
    }
    
    
    // odd generations = gametophytes (p0) just generated sporophytes
    else {
        // alternate generations.
        //old s2 - removes all individuals from pop that just reproduced (except clones)
        non_clones = p0.individuals[p0.individuals.tag != 2];
        sim.killIndividuals(non_clones);
        //reset clone tags to regular tag
        p0.individuals.tag = 1;
        
        // remove new p1s by three possible mechanisms: 
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)
        
        // 1. make a mask for random death
        if (SPO_RANDOM_DEATH_CHANCE > 0) {
            random_death = (runif(p1.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p1.individuals[random_death]);
        }
        
        // 2. maternal effect re-calculation 
        if (GAM_MATERNAL_EFFECT > 0) {
            //temp fitness scaling
            scale = SPO_POP_SIZE / p1.individualCount;
            //calculations
            maternal_fitnesses = GAM_MATERNAL_EFFECT*p1.individuals.getValue("maternal_fitness");
            child_fitnesses = (1 - GAM_MATERNAL_EFFECT)*p1.cachedFitness(NULL);
            corrected_fitnesses = scale*(maternal_fitnesses + child_fitnesses);
            
            //remove inds based on maternal effects
            to_kill = (runif(p1.individualCount) > corrected_fitnesses);
            sim.killIndividuals(p1.individuals[to_kill]);
        }
        
        
        //Set up for the next tick
        
        // fitness is scaled relative to number of inds in p1
        p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount;
    }
"""
#--------------------------------------------------------------



P0_FITNESS_SCALE_DEFAULT = "p0.fitnessScaling = GAM_POP_SIZE / p0.individualCount"
P1_FITNESS_SCALE_DEFAULT = "p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount"
WF_FITNESS_SCALE = """inds = sim.subpopulations.individuals;
    p1.fitnessScaling = K / sum(inds.fitnessScaling);"""

WF_REPRO = """
    // parents are chosen proportional to fitness
    inds = p1.individuals;
    fitness = p1.cachedFitness(NULL);
    parents1 = sample(inds, K, replace=T, weights=fitness);
    parents2 = sample(inds, K, replace=T, weights=fitness);
    for (i in seqLen(K))
        p1.addCrossed(parents1[i], parents2[i]);
    self.active = 0;
"""

EARLY_WITH_GAM_K = """
// diploids (p1) just generated haploid gametophytes
    if (community.tick % 2 == 0) {{

        // fitness affects gametophyte survival
        {p0_fitnessScaling};

       //set mutation rate for haploids
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);

        // alternate generations.
        // deactivate survival callback to remove inds from p1 subpop.
        // activate survival callback to remove inds from p0 subpop.
        s1.active = 1;
        s2.active = 0;
        
        // add random chance of death, maternal effects, or fitness calculation.
        // activate survival(p0) to set fitness method for p0s
        // deactivate survival(p1) to set fitness method for p1s
        s3.active = 1;
        s4.active = 0;
        
        // set reproduction function to be used when next generation starts.
        // deactivate reproduction(p1) diploids reproduce, haploids don't.
        // activate reproduction(p0) haploids reproduce, diploids don't.        
        s5.active = 1;
        s6.active = 0;

        {p0activate}
        {p0deactivate}
    }}


    // odd generations = gametophytes (p0) just generated sporophytes
    else {{

        // fitness affects sporophytes
        {p1_fitnessScaling};
        p0.fitnessScaling = GAM_K/ p0.individualCount;

        // set mutation rate to sporophyte rate
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);
        
        // alternate generations.
        // deactivate survival(p0) to remove inds from p0 subpop (i.e., keep p0)
        // activate survival(p1) to remove inds from p1 subpop.
        s1.active = 0;
        s2.active = 1;
        
        // add random chance of death, maternal effects, or fitness calculation.
        // deactivate survival(p0) to set fitness method for p0s
        // activate survival(p1) to set fitness method for p1s
        s3.active = 0;
        s4.active = 1;
        
        // set reproduction function to be used when next generation starts.
        // activate reproduction(p1) diploids reproduce, haploids don't.
        // deactivate reproduction(p0) haploids reproduce, diploids don't.        
        s5.active = 0;
        s6.active = 1;

        {p1activate}
        {p1deactivate}
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
    if (community.tick % 2 == 0) {{
        {muts}
    {late}
"""

SUB_MUTS = """
        fixedCount = p1.individualCount * 2 + p0.individualCount; // p1=diploid sporophytes, p0=haploid gametophytes
        mut{idx} = sim.mutationsOfType({mut});
        count{idx} = sim.mutationFrequencies(NULL, mut{idx});
        if (any(counts7 == fixedCount))
            sim.subpopulations.genomes.removeMutations(muts{idx}[counts{idx} == fixedCount], T);
"""

#-----------------------------------------------
# gam_maternal_effect
GAM_MATERNAL_EFFECT_ON_P1 = """
	// Gametophyte mother fitness affects sporophyte survival
    maternal_fitness = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_fitness)) {
        corrected_fitness = (maternal_fitness * GAM_MATERNAL_EFFECT) + fitness * (1 - GAM_MATERNAL_EFFECT);
        return (runif(1) < corrected_fitness);
    }
"""


# spo_maternal_effect
SPO_MATERNAL_EFFECT_ON_P0 = """
    // Sporophyte mother fitness affects gametophyte survival
    maternal_fitness = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_fitness)) {
        corrected_fitness = (maternal_fitness * SPO_MATERNAL_EFFECT) + fitness * (1 - SPO_MATERNAL_EFFECT);
        return (runif(1) < corrected_fitness);
    }
"""

#-----------------------------------------------
# spo_random_death_chance
# gam_random_death_chance

SURV = """
// during haploid generation, remove p1 individuals unless it is a clone.
s1 survival(p1) {{
    //clones survive to next gen and are re-tagged with parental tag
    if (individual.tag == 4) {{
            individual.tag = individual.getValue("parentid");
            return NULL;
        }}
    return F;
}}

// alternate to other generation (diploid)
// remove p0 individuals unless it is a clone.
s2 survival(p0) {{
    // allow clones to persist until next generation, tag is reset to parent tag.
    if (individual.tag == 2) {{
        individual.tag = individual.getValue("parentid");
        return NULL;
    }}
    return F;
}}

// remove p0s by three possible mechanisms: 
// 1. random chance of death.
// 2. maternal effect.
// 3. returns NULL, meaning use the standard fitness prob.
s3 survival(p0) {{
    if (runif(1) < GAM_RANDOM_DEATH_CHANCE)
        return F;
    {p1_maternal_effect}
    return NULL;
}}


// remove p1s by three possible mechanisms: 
// 1. random chance of death.
// 2. using fitness re-calculated to include maternal effect.
// 3. returns NULL, meaning use the standard fitness prob.

s4 survival(p1) {{
    if (runif(1) < SPO_RANDOM_DEATH_CHANCE)
        return F;
    {p1_maternal_effect}
    return NULL;
}}
"""

#note relFitness was replaced by `effect` in SLiM 4.0
HAP_MUT_FITNESS = """
    if (individual.subpopulation == p1)
        return 1.0;
    else
        return effect;
"""

DIP_MUT_FITNESS = """
    if (individual.subpopulation == p0)
        return 1.0;
    else
        return effect;
"""


# METADATA = """
# metadata=Dictionary(
# "spo_mutation_rate", SPO_MUTATION_RATE,
# "recombination_rate", 1e-7,
# "spo_population_size", SPO_POP_SIZE,
# "gam_population_size", GAM_POP_SIZE,
# "fixed m4 muts", sum(sim.substitutions.mutationType == m4),
# "fixed m5 muts", sum(sim.substitutions.mutationType == m5))
# """

METADATA = """
metadata=Dictionary(
"spo_mutation_rate", ifelse(SPO_MUTATION_RATE, sim.chromosome.mutationRates,
"recombination_rate", sim.chromosome.recombinationRates,
"fixed m4 muts", sum(sim.substitutions.mutationType == m4),
"fixed m5 muts", sum(sim.substitutions.mutationType == m5))
"""

DEBUG = """
// DEBUGGING
function (void)report(s$ title) {
    
    // print a title
    cat(title + "\n");
    
    // get the total number of genomes
    
    cat(format('gen=%i, ', community.tick));
    cat(format('ngenomes=%i, ', sim.subpopulations.genomesNonNull.size())); // BCH 
    cat(format('(p1=%i, ', 2 * p1.individuals.size()));
    cat(format('p0=%i, ', sum(p0.individuals.tag != 1)));
    cat(format('p0 clone=%i)\n', sum(p0.individuals.tag == 1)));
    cat("\n---------------------- fixed=" + sim.substitutions.size() + "\n\n");
    
    }
"""