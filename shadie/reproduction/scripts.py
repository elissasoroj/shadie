#!/usr/bin/env python

"""Generic SLiM4 scripts.

"""

# Standard alternation of generations early script
EARLY = """
// diploids (p1) just generated haploid gametophytes into p0
    if (community.tick % 2 == 0) {{

        // fitness affects gametophyte survival
        {p0_fitnessScaling};

        // set mutation rate for haploids
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
        // fitness is scaled relative to number of inds in p1
        {p1_fitnessScaling};

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

        {p1deactivate}
        {p1activate}
    }}
"""

# EARLY modified for the inclusion of carrying capacity in gametes
# EARLY_WITH_GAM_K = EARLY.split("\n")
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

        // fitness affects sporophytes (with haploid K capacity)
        {p1_fitnessScaling};
        p0.fitnessScaling = GAM_K / p0.individualCount;

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


# Used by EARLY_GAM_K to set {p0_fitnessScaling: P0_FITNESS_SCALE_DEFAULT}
P0_FITNESS_SCALE_DEFAULT = "p0.fitnessScaling = GAM_POP_SIZE / p0.individualCount"

# Used by EARLY_GAM_K to set {p1_fitnessScaling: P1_FITNESS_SCALE_DEFAULT}
P1_FITNESS_SCALE_DEFAULT = "p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount"

# ...
WF_FITNESS_SCALE = """inds = sim.subpopulations.individuals;
    p1.fitnessScaling = K / sum(inds.fitnessScaling);"""

# ...
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

# GAM_K
# GAM_MUTATION_RATE
# SPO_MUTATION_RATE
# SPO_POP_SIZE


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


CHECK_FOR_MUT_FIX = """
    1:{sim_time} late() {

        // only check occasionally for fixation
        if (sim.generation % 10 == 0) {

            // get the total number of genomes
            ngenomes = (2 * size(p1.individuals)) + size(p0.individuals)

            // iterate over mutationTypes in the simulation
            for (mtype in sim.mutationTypes) {
                cat('gen' + sim.generation + 'mutationType' + mtype)

                // get all mutations of this type
                muts = sim.mutationsOfType({mtype})

                //
                count = sim.mutationCounts(NULL, muts)
                freqs = sim.mutationFrequencies(NULL, muts);
                cat('count' + count)
                cat('freq' + freqs)
            }
        }
    }
"""


# this is used to check whether a substitution is completely fixed.
# A string with this chunk repeated for each mutationtype is then
# substituted into the LATE call below.
ITER_CHECK_MUT_IS_SUB = """\
        mut_{mut} = sim.mutationsOfType({mut});
        freqs = sim.mutationFrequencies(NULL, mut_{mut});

        // mutations are removed if their frequency is 0.5?
        if (any(freqs == 0.5)) {{
            sim.subpopulations.genomes.removeMutations(mut{idx}[freqs == 0.5], T);
            print('removed fixed mutation');
        }}
"""

# # this is used to check whether a substitution is completely fixed.
# ITER_CHECK_MUT_IS_SUB = """
#         count{idx} = sim.mutationFrequencies(NULL, sim.mutationsOfType({mut}));
#         if (any(counts{idx} == fixedCount))
#             sim.subpopulations.genomes.removeMutations(muts{idx}[counts{idx} == fixedCount], T);
# """

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


# DEPRECATED IN FERNS CURRENTLY FOR A MORE SPECIFIC
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

HAP_MUT_FITNESS = """
    if (individual.subpopulation == p1)
        return 1.0;
    else
        return relFitness;
"""

DIP_MUT_FITNESS = """
    if (individual.subpopulation == p0)
        return 1.0;
    else
        return relFitness;
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