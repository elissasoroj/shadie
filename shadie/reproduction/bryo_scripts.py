#!/usr/bin/env python

"""Bryophyte specific SLIM script snippets used for string substitution.

"""

# shadie DEFINITIONS
DEFS_BRYO_DIO = """
// model: dioicous bryophyte
// p0 = haploid population
// p1 = diploid population
// 1 = gametophyte (1N) tag
// 2 = gametophyte clones (1N) tag
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag; ----------- none in this model

// L0: Female = T, Male = F
"""

DEFS_BRYO_MONO = """
// model: monoicous bryophyte
// p0 = haploid population
// p1 = diploid population
// 1 = gametophyte (1N) tag
// 2 = gametophyte clones (1N) tag
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag; ----------- none in this model

// L0: Female = T, Male = F
"""

# Parameters:
# GAM_FEMALE_TO_MALE_RATIO
# -------------------------
# TAGS
# 1M-2M, 2M-3M
REPRO_BRYO_DIO_P1 = """
    ind = individual;

    // fitness-based determination of how many spores are created by this ind
    ind_fitness = p1.cachedFitness(ind.index);
    max_fitness = max(p1.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;

    // spore_vector = sample(c(0,1), SPO_SPORES_PER, replace = T, weights = c((1-ind_fitness_scaled), ind_fitness_scaled));
    // spores = sum(spore_vector);
    // sample spores weighted by fitness
    spores = rbinom(1, SPO_SPORES_PER, ind_fitness_scaled); // REPLACE ABOVE

    // each spore produces its own recombinant breakpoints
    for (rep in seqLen(spores)) {
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);

        // create four meiotic products. If later two of these mate with each other
        // it is an example of sporophytic selfing. Because we need to be able to match
        // sibling gametes at that time we tag them now with their sporophyte parent's index.
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL, parent1 = ind);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL, parent1 = ind);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL, parent1 = ind);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL, parent1 = ind);
        children = c(child1, child2, child3, child4);
        children.tag = 1; //gametophyte tag
        children.tagL0 = (runif(4) < GAM_FEMALE_TO_MALE_RATIO);

    }
"""

# Parameters:
# GAM_CLONE_RATE
# GAM_CLONES_PER
# GAM_ARCHEGONIA_PER
# SPO_SELF_RATE
# SPO_SELF_RATE_PER_EGG
# GAM_MATERNAL_EFFECT
# -------------------------
# TAGS:
# 2, 3
REPRO_BRYO_DIO_P0 = """
    // assumes archegonia only produce one egg in bryophytes
    eggs = GAM_ARCHEGONIA_PER;

    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < GAM_CLONE_RATE) {
        for (i in seqLen(GAM_CLONES_PER)) {
            child = p0.addRecombinant(individual.genome1, NULL, NULL,
            NULL, NULL, NULL, parent1 = individual); // only 1 parent recorded
            child.tag = 2; // marked as clone
            child.tagL0 = individual.tagL0; // inherits parent sex
            // gametophyte maternal effect not applicable for clones = neutral
            if (GAM_MATERNAL_EFFECT > 0)
                child.setValue("maternal_fitness", NULL);
        }
    }

    // the original plant making the clones can still reproduce regardless
    // of whether or not it cloned.
    // iterate over each egg to find a mate (self, sib, or outcross)
    // NOTE: each gametophyte gives rise to antheridia that produce thousands of
    // clonal sperm. Because of this, sperm is not removed from the mating pool when used.

    // Reproduction scripts run only female gametophytes
    if (individual.tagL0) {

        // archegonia only produce one egg in bryophytes
        eggs = GAM_ARCHEGONIA_PER;

        // In monoicous bryophytes, gametophytes are either hermaphroditic or male,
        // so "males" includes all the hermaphrodites as well as male gametophytes
        males = p0.individuals[p0.individuals.tagL0==F];

        // if selfing is possible then get all sibling males
        if (SPO_SELF_RATE_PER_EGG > 0)
            // shared parent count for sibs is > 0 (1 or 2)
            siblings = males[individual.sharedParentCount(males)!=0];

        for (rep in 1:eggs) {

            // idea: random_chance_of ... here?
            // weighted sampling: each egg gets gam-selfed, spo-selfed, or outcrossed.
            // Note: shadie enforces that GAM_SELF_RATE + SPO_SELF_RATE !> 1
            mode = sample(
                x=c(1, 2, 3),
                size=1,
                weights=c(0, SPO_SELF_RATE_PER_EGG, 1 - (SPO_SELF_RATE_PER_EGG))
                );

            // intra-gametophytic selfed
            if (mode == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL,
                individual.genome1, NULL, NULL, parent1 = individual, parent2 = individual);
                child.tag = 3; //sporophyte tag
                // gametophyte maternal effect on new sporophyte
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }

            // inter-gametophytic selfed individual (same sporo parent)
            // only occurs IF a sibling gametophyte is still alive.
            else if (mode == 2) {
                if (siblings.size() > 0) {
                    sibling = sample(siblings, 1);
                    child = p1.addRecombinant(individual.genome1, NULL, NULL,
                    sibling.genome1, NULL, NULL, parent1 = individual, parent2=sibling);
                    child.tag = 3; //sporophyte tag
                    // gametophyte maternal effect on new sporophyte
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
            // outcrossing individual samples any other p0, and checks that it is outcrossing
            // only occurs if a non-sib gametophyte is still alive.

            else {
                outcross_sperms = males[individual.sharedParentCount(males)==0];
                if (! isNULL(outcross_sperms)) {
                    sperm = sample(outcross_sperms, 1);
                    child = p1.addRecombinant(individual.genome1, NULL, NULL,
                    sperm.genome1, NULL, NULL, parent1 = individual, parent2=sperm);
                    child.tag = 3; //sporophyte tag
                    // gametophyte maternal effect on new sporophyte
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
    }
"""


# ==PARAMETERS==
# SPO_SPORES_PER
# GAM_FEMALE_TO_MALE_RATIO
# ------------------
# TAGS
REPRO_BRYO_MONO_P1 = """
    ind = individual;

    // fitness-based determination of how many spores are created by this ind
    ind_fitness = p1.cachedFitness(ind.index);
    max_fitness = max(p1.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;

    //spore_vector = sample(c(0,1), SPO_SPORES_PER, replace = T, weights = c((1-ind_fitness_scaled), ind_fitness_scaled));
    //spores = sum(spore_vector);
    spores = rbinom(1, SPO_SPORES_PER, ind_fitness_scaled); // REPLACE ABOVE

    // each spore produces its own recombinant breakpoints
    for (rep in seqLen(spores)) {
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);

        // create four meiotic products. If later two of these mate with each other
        // it is an example of sporophytic selfing. Because we need to be able to match
        // sibling gametes at that time we tag them now with their sporophyte parent's index.
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL, parent1 = ind);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL, parent1 = ind);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL, parent1 = ind);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL, parent1 = ind);
        children = c(child1, child2, child3, child4);
        children.tag = 1; //gametophyte tag
        children.tagL0 = (runif(4) < GAM_FEMALE_TO_MALE_RATIO);
    }
"""

# ==PARAMETERS==
# GAM_ARCHEGONIA_PER
# GAM_CLONE_RATE
# GAM_MATERNAL_EFFECT
# ------------------
# TAGS
# 2, >1M-
REPRO_BRYO_MONO_P0 = """
    // assumes archegonia only produce one egg in bryophytes
    eggs = GAM_ARCHEGONIA_PER;

    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < GAM_CLONE_RATE) {
        for (i in seqLen(GAM_CLONES_PER)) {

            // only 1 parent recorded
            child = p0.addRecombinant(
                individual.genome1, NULL, NULL, NULL, NULL, NULL, parent1=individual);

            // marked as clone
            child.tag = 2;

            // inherits parent sex
            child.tagL0 = individual.tagL0;

            // gametophyte maternal effect not applicable for clones = neutral
            if (GAM_MATERNAL_EFFECT > 0)
                child.setValue("maternal_fitness", NULL);
        }
    }

    // the original plant making the clones can still reproduce regardless
    // of whether or not it cloned.
    // iterate over each egg to find a mate (self, sib, or outcross)
    // NOTE: each gametophyte gives rise to antheridia that produce thousands of
    // clonal sperm. Because of this, sperm is not removed from the mating pool when used.

    // Reproduction scripts run only female gametophytes
    if (individual.tagL0) {

        // archegonia only produce one egg in bryophytes
        eggs = GAM_ARCHEGONIA_PER;

        // In monoicous bryophytes, gametophytes are either hermaphroditic or male,
        // so "males" includes all the hermaphrodites as well as male gametophytes
        males = p0.individuals;

        // if selfing is possible then get all sibling males
        if (SPO_SELF_RATE_PER_EGG > 0)
            // shared parent count for sibs is > 0 (1 or 2)
            siblings = males[individual.sharedParentCount(males)!=0];

        for (rep in 1:eggs) {

            // idea: random_chance_of ... here?
            // weighted sampling: each egg gets gam-selfed, spo-selfed, or outcrossed.
            // Note: shadie enforces that GAM_SELF_RATE + SPO_SELF_RATE !> 1
            mode = sample(
                x=c(1, 2, 3),
                size=1,
                weights=c(GAM_SELF_RATE_PER_EGG, SPO_SELF_RATE_PER_EGG, 1 - (GAM_SELF_RATE_PER_EGG + SPO_SELF_RATE_PER_EGG))
            );

            // intra-gametophytic selfed
            if (mode == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL,
                individual.genome1, NULL, NULL, parent1 = individual, parent2 = individual);
                child.tag = 3; //sporophyte tag
                //gametophyte maternal effect on new sporophyte
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }

            // inter-gametophytic selfed individual (same sporo parent)
            // only occurs IF a sibling gametophyte is still alive.
            else if (mode == 2) {
                if (siblings.size() > 0) {
                    sibling = sample(siblings, 1);
                    child = p1.addRecombinant(individual.genome1, NULL, NULL,
                    sibling.genome1, NULL, NULL, parent1 = individual, parent2 = sibling);
                    child.tag = 3; //sporophyte tag
                    //gametophyte maternal effect on new sporophyte
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
            // outcrossing individual samples any other p0, and checks that it is outcrossing
            // only occurs if a non-sib gametophyte is still alive.

            else {
                outcross_sperms = males[individual.sharedParentCount(males)==0];
                if (! isNULL(outcross_sperms)) {
                    sperm = sample(outcross_sperms, 1);
                    child = p1.addRecombinant(individual.genome1, NULL, NULL,
                    sperm.genome1, NULL, NULL, parent1 = individual, parent2=sperm);
                    child.tag = 3; //sporophyte tag
                    //gametophyte maternal effect on new sporophyte
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
    }
"""
