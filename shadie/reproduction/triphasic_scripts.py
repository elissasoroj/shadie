#!/usr/bin/env python

"""Generic SLiM4 scripts for triphasic models.

"""

######################################################################
# standard first script

# ==PARAMETERS==
# SPO_MUTATION_RATE
# GAM_MUTATION_RATE
###################
FIRST = """
    if (sim.cycle % 3 == 2) {
        // reproduction(p2) will create GAMETOPHYTES into p1.

        // set reproduction function to be used next generation
        s0.active = 0; //CSPO
        s1.active = 0; // GAM
        s2.active = 1; // TSPO

        // set mutation rate of offspring
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);
    }

    else if (sim.cycle % 3 == 1) {
        // reproduction(p1) will create TETRASPOROPHYTES in p2.

        // set reproduction function to be used next generation
        s0.active = 1; //CSPO
        s1.active = 0; // GAM
        s2.active = 0; // TSPO

        // set mutation rate of offspring
        sim.chromosome.setMutationRate(TSPO_MUTATION_RATE);
    }

    else { // %3 == 0
        // reproduction(p1) will create CARPOSPOROPHYTES in p0.

        // set reproduction function to be used next generation
        s0.active = 0; //CSPO
        s1.active = 1; // GAM
        s2.active = 0; // TSPO

        // set mutation rate of offspring
        sim.chromosome.setMutationRate(CSPO_MUTATION_RATE);
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
    if (community.tick % 3 == 2) {{
        // reproduction(p2) just finished creating GAMETOPHYTES into p1.
        // this generation will calculate fitness and apply selection on them.
        
        {tetrasporophyte_clones}

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

        // 2. No tetrasporophyte maternal effect on p1

        // 3. Set up for next tick
        {p1_survival_effects}
    }}


    else if (community.tick % 3 == 1) {{
        // reproduction(p0) just finished creating TETRASPOROPHYTES into p2.
        // this generation will calculate fitness and apply selection on them.

        sim.killIndividuals(p0.individuals);

        // new p1s removed by three possible mechanisms:
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)

        // 1. make a mask for random death
        if (GAM_RANDOM_DEATH_CHANCE > 0) {{
            random_death = (runif(p1.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p1.individuals[random_death]);
        }}

        // 2. No carposporophyte maternal effect on p2

        // 3. Set up for next tick
        {p2_survival_effects}

    }}

    else {{
        //reproduction (p1) just finished creating CARPOSPOROPHYTES into p0
        // this generation will calculate fitness and apply selection on them.
        
        // alternate generations.
            {gametophyte_clones}

        // remove new p2s by three possible mechanisms:
        // 1. random chance of death.
        // 2. using fitness re-calculated to include maternal effect.
        // 3. individual goes on next tick (unless fitness kills it)

        // 1. make a mask for random death
        if (CSPO_RANDOM_DEATH_CHANCE > 0) {{
            random_death = (runif(p2.individualCount) < GAM_RANDOM_DEATH_CHANCE);
            sim.killIndividuals(p2.individuals[random_death]);
        }}

        {gam_maternal_effect}

        // 3. Set up for next tick
        //all carposporophytes survive; populations is controlled
        //by gametophyte population and GAM_CARPOGONIA_PER param
    }}
"""


###############################################################
# cloning scripts

TSPO_CLONES = """
    // remove parents, except for clones
    non_clones = p2.individuals[p2.individuals.tag != 4];
    sim.killIndividuals(non_clones);
    // reset clone tags to regular tag
    p2.individuals.tag = 3;
"""

NO_TSPO_CLONES = """
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
P2_FITNESS_SCALE_DEFAULT = "p2.fitnessScaling = TSPO_POP_SIZE / p2.individualCount;"

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

FITNESS_AFFECTS_TSPO_REPRODUCTION = """
// fitness-based determination of how many tetraspores are created by this ind
    ind_fitness = p2.cachedFitness(ind.index);
    max_fitness = max(p2.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;

    tetraspores = rbinom(1, TSPO_MAX_SPORES_PER, ind_fitness_scaled);
"""

CONSTANT_TSPORES = "tetraspores = TSPO_MAX_SPORES_PER;"


FITNESS_AFFECTS_SPERM_SUCCESS = """
// fitness-based determination of sperm sampling
    sperm_fitness_vector = p1.cachedFitness(outcross_sperms.index);
    sperm = sample(outcross_sperms, 1, weights = sperm_fitness_vector);
"""

RANDOM_MATING = "sperm = sample(outcross_sperms, 1);"


#################################################################
# gam_maternal_effect
GAM_MATERNAL_EFFECT_ON_P0 = """
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
