#!/usr/bin/env python

"""
Pteridophyte specific SLIM script snippets used for string substitution.
"""
PTER_FITNESS_SCALE = """p0.fitnessScaling = (GAM_POP_SIZE * 
    (1 / GAM_FEMALE_TO_MALE_RATIO)) * (2 / GAM_ARCHEGONIA_PER) / p0.individualCount;"""

# PARAMETERS
# -------------------------
# TAGS
# 2, 3, 4, >1M

DEFS_PTER_HOMOSPORE = """
// model: homosporous pteridophyte
// p0 = haploid population
// p1 = diploid population
// >1M- tmp sex tags
// <2000000 = hermaphroditic gametophyte (1N)
// >2000000 = male gametophyte (1N)
// 2 = gametophyte clones (1N) tag
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag
"""
# ==PARAMETERS==
#SPO_CLONE_RATE
#SPO_CLONES_PER
#SPO_SPORES_PER
#GAM_FEMALE_TO_MALE_RATIO
# -------------------------
# TAGS
# 4, >1M-
REPRO_PTER_HOMOSPORE_P1 = """
    ind = individual;
    
    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < SPO_CLONE_RATE) {
        for (i in 1:SPO_CLONES_PER) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome2, NULL, NULL);
            child.tag = 4;
            child.setValue("parentid", individual.tag);
        }
    }
    
    // parent tag is 1M + the parent's unique index
    //hermaphrodite tags start with 1
    htag = 1000000 + ind.index;
    //male tags start with 2
    mtag = 2000000 + ind.index;
    
    // fitness-based determination of how many spores are created by this ind
    ind_fitness = p1.cachedFitness(ind.index);
    max_fitness = max(p1.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness; 
    
    spore_vector = sample(c(0,1), SPO_SPORES_PER, replace = T, weights = c((1-ind_fitness_scaled), ind_fitness_scaled));
    spores = sum(spore_vector);

    // each spore produces its own recombinant breakpoints 
    for (rep in 1:spores) {
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        
        // create four meiotic products. If later two of these mate with each other
        // it is an example of sporophytic selfing. Because we need to be able to match
        // sibling gametes at that time we tag them now with their sporophyte parent's index.
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
        child1.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, htag, mtag);
        child2.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, htag, mtag);
        child3.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, htag, mtag);
        child4.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, htag, mtag);
        children = c(child1, child2, child3, child4);
        if (SPO_MATERNAL_EFFECT > 0)
            children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
    }
"""

# PARAMETERS
#GAM_CLONE_RATE
#GAM_CLONES_PER
#GAM_MATERNAL_EFFECT
#SPO_SELF_RATE_PER_EGG
#SPO_SELF_RATE
# -------------------------
# TAGS
# 2, >1M-
REPRO_PTER_HOMOSPORE_P0 = """
    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < GAM_CLONE_RATE) {
        for (i in 1:GAM_CLONES_PER) {
            child = p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL);
            child.tag = 2;
            child.setValue("parentid", individual.tag);
            if (GAM_MATERNAL_EFFECT > 0)
                child.setValue("maternal_fitness", NULL);
        }
    }
    
    // iterate over each egg to find a mate (self, sib, or outcross)
    // NOTE: each gametophyte gives rise to antheridia that produce thousands of
    // clonal sperm. Because of this, sperm is not removed from the mating pool when used. 
    
    //Reproduction scripts run only hermaphroditic gametophytes
    if (individual.tag < 2000000) {
        
        // get all males that could fertilize an egg of this hermaphrodite
        //In homosporous ferns, most gametophytes are hermaphroditic or male, 
        //so "males" includes all the hermaphrodites as well as male gametophytes
        males = p0.individuals; 
        
        // if selfing is possible then get all sibling males
        if (SPO_SELF_RATE_PER_EGG > 0)
            siblings = males[males.tag == 1000000 + individual.tag];
        
        // iterate over each reproductive opportunity (archegonia) in this hermaphrodite.
        // A hermaphrodite could produce multiple eggs per arch, but they are identical, 
        // so we do not bother to model that for now ...
        for (rep in 1:GAM_ARCHEGONIA_PER) {
            
            // weighted sampling: each egg gets gam-selfed, spo-selfed, or outcrossed.
            // Note: shadie enforces that GAM_SELF_RATE + SPO_SELF_RATE !> 1
            mode = sample(
                x=c(1, 2, 3),
                size=1,
                weights=c(GAM_SELF_RATE_PER_EGG, SPO_SELF_RATE_PER_EGG, 1 - (GAM_SELF_RATE_PER_EGG + SPO_SELF_RATE_PER_EGG))
                );
            
            // intra-gametophytic selfed
            if (mode == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
                child.tag = 3;
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
            
            // inter-gametophytic selfed individual (same sporo parent)
            // only occurs IF a sibling gametophyte is still alive.
            else if (mode == 2) {
                // sibling = p0.sampleIndividuals(1, tag=1000000 + individual.tag);
                if (siblings.size() > 0) {
                    sibling = sample(siblings, 1);
                    child = p1.addRecombinant(individual.genome1, NULL, NULL, sibling.genome1, NULL, NULL);
                    child.tag = 3;
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
            
            // outcrossing individual samples any p0 that is not same tag.
            // only occurs if a non-sib gametophyte is still alive.
            else {
                // try at most 10 times to find a non-sib sperm, then skip.
                // males = p0.individuals[p0.individuals.tag > 2000000];
                for (trial in 1:10) {
                    sperm = sample(males, 1);
                    if (sperm.tag != individual.tag + 1000000) {
                        child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                        child.tag = 3;
                        if (GAM_MATERNAL_EFFECT > 0)
                            child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                        break;
                    }
                }
            }
        }
    }
"""

# PARAMETERS
# -------------------------
# TAGS
# 2, 3, 4, >1M-
#-----------------------------------------------------------------------
DEFS_PTER_HETEROSPORE = """
// model: heterosporous pteridophyte
// p0 = haploid population
// p1 = diploid population
// >1M- tmp sex tags
// <2000000 = female gametophyte (1N)
// >2000000 = male gametophyte (1N)
// 2 = gametophyte clones (1N) tag -------------- None in this model
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag
"""


# ==PARAMETERS==
#SPO_CLONE_RATE
#SPO_CLONES_PER
#GAM_FEMALE_TO_MALE_RATIO
#SPO_MATERNAL_EFFECT
#SPO_SPORES_PER
# -------------------------
# TAGS
# 4, >1M-
REPRO_PTER_HETEROSPORE_P1 = """
     ind = individual;
    
    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < SPO_CLONE_RATE) {
        for (i in 1:SPO_CLONES_PER) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome2, NULL, NULL);
            child.tag = 4;
            child.setValue("parentid", individual.tag);
        }
    }
    
    // parent tag is 1M + the parents unique index
    //female tags start with 1
    ftag = 1000000 + ind.index;
    //male tags start with 2
    mtag = 2000000 + ind.index;
    
    // each spore produces its own recombinant breakpoints 
    for (rep in 1:SPO_SPORES_PER) {
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        
        // create four meiotic products. If later two of these mate with each other
        // it is an example of sporophytic selfing. Because we need to be able to match
        // sibling gametes at that time we tag them now with their sporophyte parent's index.
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
        child1.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        child2.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        child3.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        child4.tag = ifelse(runif(1) < GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        children = c(child1, child2, child3, child4);
        if (SPO_MATERNAL_EFFECT > 0)
            children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
    }
"""

# PARAMETERS
#v
#GAM_ARCHEGONIA_PER
# -------------------------
# TAGS
#3, >1M-
REPRO_PTER_HETEROSPORE_P0 = """
    // iterate over each egg to find a mate (self, sib, or outcross)
    // NOTE: each gametophyte gives rise to antheridia that produce thousands of
    // clonal sperm. Because of this, sperm is not removed from the mating pool when used. 
    
    //Reproduction scripts run only female gametophytes
    if (individual.tag < 2000000) {
        
        // get all males that could fertilize an egg of this female
        males = p0.individuals[p0.individuals.tag > 2000000];
        
        // if selfing is possible then get all sibling males
        if (SPO_SELF_RATE_PER_EGG > 0)
            siblings = males[males.tag == 1000000 + individual.tag];
        
        // iterate over each reproductive opportunity (archegonia) in this female.
        // A female could produce multiple eggs per arch, but they are identical, 
        // so we do not bother to model that for now ...
        for (rep in 1:GAM_ARCHEGONIA_PER) {
            
            
            // inter-gametophytic selfed individual (same sporo parent)
            // only occurs IF a sibling gametophyte is still alive.
            if (runif(1) < SPO_SELF_RATE_PER_EGG) {
                // sibling = p0.sampleIndividuals(1, tag=1000000 + individual.tag);
                if (siblings.size() > 0) {
                    sibling = sample(siblings, 1);
                    child = p1.addRecombinant(individual.genome1, NULL, NULL, sibling.genome1, NULL, NULL);
                    child.tag = 3;
                }
            }
            
            // outcrossing individual samples any p0 that is not same tag.
            // only occurs if a non-sib gametophyte is still alive.
            else {
                // try at most 10 times to find a non-sib sperm, then skip.
                // males = p0.individuals[p0.individuals.tag > 2000000];
                for (trial in 1:10) {
                    sperm = sample(males, 1);
                    if (sperm.tag != individual.tag + 1000000) {
                        child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                        child.tag = 3;
                        break;
                    }
                }
            }
        }
    }
"""