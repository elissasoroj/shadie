#!/usr/bin/env python

"""Generic SLiM4 scripts.

"""

######################################################################
# standard first script

# ==PARAMETERS==
# SPO_MUTATION_RATE
# GAM_MUTATION_RATE
###################
FIRST = """
// alternate generations
    if (sim.cycle % 2 == 1) {
        // reproduction(p1) will create SPOROPHYTES in p2.

        // set reproduction function to be used this generation
        s1.active = 1; // GAM
        s2.active = 0; // SPO

        // set mutation rate that will apply to the offspring
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);
    }

    else {
        // reproduction(p2) will produce GAMETOPHYTES into p1.

        // set reproduction function to be used this generation
        s1.active = 0;  // GAM
        s2.active = 1;  // SPO

        // set mutation rate that will apply to the offspring
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);
    }
"""

#####################################################################
# standard early script

# ==PARAMETERS==
# GAM_RANDOM_DEATH_CHANCE
# GAM_CEILING
# GAM_POP_SIZE
# SPO_RANDOM_DEATH_CHANCE
# SPO_POP_SIZE
# -------------------------
EARLY = """
// diploids (p2) just generated haploid gametophytes into p1
    if (community.tick % 2 == 0) {{
        
        // alternate generations.
            {sporophyte_clones}

        // new p1s removed by three possible mechanisms:
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)

        // 1. make a mask for random death
        if (GAM_RANDOM_DEATH_CHANCE > 0) {{
            random_death = (runif(p1.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p1.individuals[random_death]);
        }}

        // random death also occurs to implement GAM_CEILING
        if (p1.individualCount > GAM_CEILING) {{
            to_kill = GAM_CEILING - GAM_POP_SIZE;
            death_chance = to_kill/p1.individualCount;
            random_death = sample(c(F,T), p1.individualCount, T, c(1-death_chance, death_chance));
            sim.killIndividuals(p1.individuals[random_death]);
            }}

        {spo_maternal_effect}

        // 3. Set up for next tick
        // fitness affects gametophyte survival
        {p1_survival_effects}
    }}


    // odd generations = gametophytes (p1) just generated sporophytes
    else {{
        // alternate generations.
            {gametophyte_clones}

        // remove new p2s by three possible mechanisms:
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)

        // 1. make a mask for random death
        if (SPO_RANDOM_DEATH_CHANCE > 0) {{
            random_death = (runif(p2.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p2.individuals[random_death]);
        }}
        {gam_maternal_effect}

        // 3. Set up for next tick
        // fitness is scaled relative to number of inds in p2
        {p2_survival_effects}
    }}
"""


###############################################################
# cloning scripts

SPO_CLONES = """
    // remove parents, except for clones
    non_clones = p2.individuals[p2.individuals.tag != 4];
    sim.killIndividuals(non_clones);
    // reset clone tags to regular tag
    p2.individuals.tag = 3;
"""

NO_SPO_CLONES = """
    // remove parents
    sim.killIndividuals(p2.individuals);
"""

GAM_CLONES = """
    // removes parents, except clones
    non_clones = p1.individuals[p1.individuals.tag != 2];
    sim.killIndividuals(non_clones);
    // reset clone tags to regular tag
    p1.individuals.tag = 1;
"""

NO_GAM_CLONES = """
    // remove parents
    sim.killIndividuals(p1.individuals);
"""

###############################################################
# fitness scripts
#
# ==PARAMETERS==
# GAM_POP_SIZE
# SPO_POP_SIZE
# ----------------------
P1_FITNESS_SCALE_DEFAULT = "p1.fitnessScaling = GAM_POP_SIZE / p1.individualCount;"
P2_FITNESS_SCALE_DEFAULT = "p2.fitnessScaling = SPO_POP_SIZE / p2.individualCount;"
WF_FITNESS_SCALE = """inds = sim.subpopulations.individuals;
    p2.fitnessScaling = K / sum(inds.fitnessScaling);"""

P1_RANDOM_SURVIVAL = """
    //exact population is maintained
    num_to_kill = p1.individualCount - GAM_POP_SIZE;
    if (num_to_kill > 0) {
        random_inds = sample(p1.individuals, num_to_kill);
        sim.killIndividuals(p1.individuals[random_inds.index]);           
    }
"""
P2_RANDOM_SURVIVAL = """
    //exact population is maintained
    num_to_kill = p2.individualCount - SPO_POP_SIZE;
    if (num_to_kill > 0) {
        random_inds = sample(p2.individuals, num_to_kill);
        sim.killIndividuals(p2.individuals[random_inds.index]);           
    }
"""

# -----------------------
# ==SCRIPTS FOR SOFT SELECTION==

FITNESS_AFFECTS_SPO_REPRODUCTION = """
// fitness-based determination of how many spores are created by this ind
    ind_fitness = p2.cachedFitness(ind.index);
    max_fitness = max(p2.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;

    spores = rbinom(1, SPO_SPORES_PER, ind_fitness_scaled);
"""

FITNESS_AFFECTS_GAM_MATING = """
// fitness-based determination of sperm sampling
    sperm_fitness_vector = p1.cachedFitness(outcross_sperms.index);
    sperm = sample(outcross_sperms, 1, weights = sperm_fitness_vector);
"""

WF_REPRO_SOFT = """
    // parents are chosen proportional to fitness
    inds = p2.individuals;
    fitness = p2.cachedFitness(NULL);
    parents1 = sample(inds, K, replace=T, weights=fitness);
    parents2 = sample(inds, K, replace=T, weights=fitness);
    for (i in seqLen(K))
        p2.addCrossed(parents1[i], parents2[i]);
    self.active = 0;
"""

# ==SCRIPTS FOR HARD SELECTION==

RANDOM_MATING = "sperm = sample(outcross_sperms, 1);"

CONSTANT_SPORES = "spores = SPO_SPORES_PER;"

WF_REPRO_HARD = """
    //parents chosen at random for mating
    // child number is poisson draw w/ lambda=2
    num_offspring = rpois(N, {lambda_pois}); 

    //adjust parents to maintain K
    if (p2.individualCount < N){{
        add_num = N-p2.individualCount;
        add_par = sample(p2.individuals, add_num);
        parents = c(p2.individuals, add_par);
    }}
    else
        parents = p2.individuals;

    //create the gamete pool based on poisson sampling
    gamete_pool = repEach(parents, num_offspring);

    //randomize the gamete pool
    rand_gamete_pool = sample(gamete_pool, length(gamete_pool), replace=F);
    
    //split the gametes into two parent vectors
    matings = asInteger(floor(length(gamete_pool)/2));
    parents1 = gamete_pool[0:matings];
    end = length(gamete_pool) - 1;
    parents2 = gamete_pool[matings:end];
    
    //make the offspring
    for (i in seqLen(matings))
       p2.addCrossed(parents1[i], parents2[i]);
    self.active = 0;
"""

#implemented directly in base
WF_EARLY_HARD = """
//kill parents first
    sim.killIndividuals(p2.individuals[p2.individuals.age > 0]);
"""

#goes in late() call
WF_SELECTION = """
    //enforces constant population size
    numtokill = length(p2.individuals) - N;
    if (numtokill > 0){
        fitness = p2.cachedFitness(NULL);
        fitness_weights = (1+max(fitness))-fitness;
        tokill = sample(p2.individuals, numtokill, weights = fitness_weights);
        sim.killIndividuals(tokill);
    }
"""

#################################################################
# activate/deactivate fitness callbacks

ACTIVATE = "{idx}.active = 1;"

DEACTIVATE = "{idx}.active = 0;"

#################################################################
# gam_maternal_effect
GAM_MATERNAL_EFFECT_ON_P2 = """
    // 2. maternal effect re-calculation
    if (GAM_MATERNAL_EFFECT > 0) {
        // temp fitness scaling
        scale = SPO_POP_SIZE / p2.individualCount;

        // calculations
        maternal_fitnesses = GAM_MATERNAL_EFFECT*p2.individuals.getValue("maternal_fitness");
        child_fitnesses = (1 - GAM_MATERNAL_EFFECT)*p2.cachedFitness(NULL);
        corrected_fitnesses = scale*(maternal_fitnesses + child_fitnesses);

        // remove inds based on maternal effects
        to_kill = (runif(p2.individualCount) > corrected_fitnesses);
        sim.killIndividuals(p2.individuals[to_kill]);
    }
"""

NO_GAM_MATERNAL_EFFECT = """
    // 2. No gametophyte maternal effect on p2
"""

# spo_maternal_effect
SPO_MATERNAL_EFFECT_ON_P1 = """
    // 2. maternal effect re-calculation
    if (SPO_MATERNAL_EFFECT > 0) {
        // temp fitness scaling
        scale = SPO_POP_SIZE / p2.individualCount;
        // calculations
        maternal_fitnesses = SPO_MATERNAL_EFFECT*p1.individuals.getValue("maternal_fitness");
        child_fitnesses = (1 - SPO_MATERNAL_EFFECT)*p1.cachedFitness(NULL);
        corrected_fitnesses = scale * (maternal_fitnesses + child_fitnesses);

        // remove inds based on maternal effects
        to_kill = (runif(p1.individualCount) > corrected_fitnesses);
        sim.killIndividuals(p1.individuals[to_kill]);
    }
"""

NO_SPO_MATERNAL_EFFECT = """
        // 2. No sporophyte maternal effect on p1
"""

##################################################################
# note relFitness was replaced by `effect` in SLiM 4.0
HAP_MUT_FITNESS = """
    if (individual.subpopulation == p2)
        return 1.0;
    else
        return effect;
"""

DIP_MUT_FITNESS = """
    if (individual.subpopulation == p1)
        return 1.0;
    else
        return effect;
"""
##################################################################

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
    cat(format('(p2=%i, ', 2 * p2.individuals.size()));
    cat(format('p1=%i, ', sum(p1.individuals.tag != 1)));
    cat(format('p1 clone=%i)\n', sum(p1.individuals.tag == 1)));
    cat("\n---------------------- fixed=" + sim.substitutions.size() + "\n\n");

    }
"""
