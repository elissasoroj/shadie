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
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);
    }
    
    else {  //(sim.cycle % 2 == 0)
        // reproduction(p1) will produce GAMETOPHYTES into p0.
        
        // set reproduction function to be used this generation  
        s0.active = 0; //GAM
        s1.active = 1; //SPO
        
        // set mutation rate that will apply to the offspring
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);
    }
"""
#--------------------------------------------------------------
#standard early script

# ==PARAMETERS==
# GAM_RANDOM_DEATH_CHANCE
# GAM_CEILING
# GAM_POP_SIZE
# SPO_RANDOM_DEATH_CHANCE
# SPO_POP_SIZE
# -------------------------
EARLY = """
// diploids (p1) just generated haploid gametophytes into p0
    if (community.tick % 2 == 0) {{
        // alternate generations.
        
        {sporophyte_clones}
        
        // new p0s removed by three possible mechanisms: 
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)
        
        // 1. make a mask for random death
        if (GAM_RANDOM_DEATH_CHANCE > 0) {{
            random_death = (runif(p0.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p0.individuals[random_death]);
        }}

        //random death also occurs to implement GAM_CEILING
        if (p0.individualCount > GAM_CEILING) {{
            to_kill = GAM_CEILING - GAM_POP_SIZE;
            death_chance = to_kill/p0.individualCount;
            random_death = sample(c(F,T), p0.individualCount, T, c(1-death_chance, death_chance));
            sim.killIndividuals(p0.individuals[random_death]);
            }}
        
       {spo_maternal_effect}
        
        //3. Set up for next tick
        // fitness affects gametophyte survival
        {p0_fitnessScaling}
    }}
    
    
    // odd generations = gametophytes (p0) just generated sporophytes
    else {{
        // alternate generations.
        {gametophyte_clones}
        
        // remove new p1s by three possible mechanisms: 
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)
        
        // 1. make a mask for random death
        if (SPO_RANDOM_DEATH_CHANCE > 0) {{
            random_death = (runif(p1.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p1.individuals[random_death]);
        }}
        
        {gam_maternal_effect}
        
        
        //3. Set up for next tick
        // fitness is scaled relative to number of inds in p1
        {p1_fitnessScaling}
    }}
"""
#--------------------------------------------------------------

SPO_CLONES = """
    // remove parents, except for clones
    non_clones = p1.individuals[p1.individuals.tag != 4];
    sim.killIndividuals(non_clones);
    //reset clone tags to regular tag
    p1.individuals.tag = 3;
"""

NO_SPO_CLONES = """
    // remove parents
    sim.killIndividuals(p1.individuals);
"""

GAM_CLONES = """
// removes parents, except clones
        non_clones = p0.individuals[p0.individuals.tag != 2];
        sim.killIndividuals(non_clones);
        //reset clone tags to regular tag
        p0.individuals.tag = 1;
"""

NO_GAM_CLONES = """
    // remove parents
    sim.killIndividuals(p0.individuals);
"""

#--------------------------------------------------------------

# ==PARAMETERS==
# GAM_POP_SIZE
# SPO_POP_SIZE
# ----------------------
P0_FITNESS_SCALE_DEFAULT = "p0.fitnessScaling = GAM_POP_SIZE / p0.individualCount;"
P1_FITNESS_SCALE_DEFAULT = "p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount;"
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
#-----------------------------------------------
#activate/deactivate fitness callbacks

ACTIVATE = "{idx}.active = 1;"

DEACTIVATE = "{idx}.active = 0;"

#-----------------------------------------------
# gam_maternal_effect
GAM_MATERNAL_EFFECT_ON_P1 = """
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
"""

NO_GAM_MATERNAL_EFFECT = """
    // 2. No gametophyte maternal effect on P1
"""

# spo_maternal_effect
SPO_MATERNAL_EFFECT_ON_P0 = """
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
"""

NO_SPO_MATERNAL_EFFECT = """
    // 2. No sporophyte maternal effect on P0
"""

#-----------------------------------------------
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
#-----------------------------------------------

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