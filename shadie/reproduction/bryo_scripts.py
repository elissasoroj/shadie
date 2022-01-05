#!/usr/bin/env python

"""
Bryophyte specific SLIM script snippets used for string substitution.
"""
DEFS_BRYO = """
// shadie DEFINITIONS
// p0 = haploid population
// p1 = diploid population
// >2000000 = male gametophyte (1N)tag
// <2000000 = female gametophyte (1N)tag
// 2 = gametophyte clones (1N)tag
// 3 = sporophyte (2N)tag
// >1M- tmp tags
// Note: sporophytes cannot clone here.
"""

# Parameters:
#GAM_FEMALE_TO_MALE_RATIO
# -------------------------
# TAGS
# 1M-2M, 2M-3M
REPRO_BRYO_DIO_P1 = """
    ind = individual;
    
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
        p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL).tag=ifelse(runif(1)<GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL).tag=ifelse(runif(1)<GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL).tag=ifelse(runif(1)<GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL).tag=ifelse(runif(1)<GAM_FEMALE_TO_MALE_RATIO, ftag, mtag);
        }
"""

# Parameters:
# GAM_CLONE_RATE
# GAM_CLONES_PER
# GAM_ARCHEGONIA_PER
# GAM_SELF_RATE
# GAM_SELF_RATE_PER_EGG
# SPO_SELF_RATE
# SPO_SELF_RATE_PER_EGG
# GAM_MATERNAL_EFFECT
# -------------------------
# TAGS:
# 2, 3
REPRO_BRYO_DIO_P0 = """
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
    
    // the original plant making the clones can still reproduce regardless 
    // of whether or not it cloned.
    // iterate over each egg to find a mate (self, sib, or outcross)
    // NOTE: each gametophyte gives rise to antheridia that produce thousands of
    // clonal sperm. Because of this, sperm is not removed from the mating pool when used. 
    
    //reproduction scripts run only on female gametophytes
    if (individual.tag < 2000000) {
        
        // archegonia only produce one egg in bryophytes
        eggs = GAM_ARCHEGONIA_PER;
        
        for (rep in 1:eggs) {
            
            // idea: random_chance_of ... here?
            // weighted sampling: each egg gets gam-selfed, spo-selfed, or outcrossed.
            // Note: shadie enforces that GAM_SELF_RATE + SPO_SELF_RATE !> 1
            mode = sample(
                x=c(1, 2, 3),
                size=1,
                weights=c(0.0, SPO_SELF_RATE_PER_EGG, 1 - (GAM_SELF_RATE_PER_EGG + SPO_SELF_RATE_PER_EGG))
                );
            
            // intra-gametic selfed (only possible in monoicous)
            if (mode == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
                child.tag = 3;
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
            
            // inter-gametophytic selfed individual (same sporo parent)
            // only occurs IF a sibling gametophyte is still alive.
            else if (mode == 2) {
                sibling = p0.sampleIndividuals(1, tag=1000000+individual.tag);
                if (size(sibling)==1) {
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
                males = p0.individuals[p0.individuals.tag > 2000000];
                for (trial in 1:10) {
                    sperm = sample(males, 1);
                    if (sperm.tag != individual.tag+1000000) {
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


## PARAMETERS
# spo_spores_per
# ------------------
# TAGS
# 0, 4, 5
REPRO_BRYO_MONO_P1 = """
    if (individual.tag == 0){
        meiosis_reps = asInteger(SPO_SPORES_PER/4);
        make_spores(individual, meiosis_reps);
    }
    
    //sporophytic selfing
    if (individual.tag == 5) //sporophyte selfs
        sporophyte_selfs(individual);
"""

## PARAMETERS
# 
# ------------------
# TAGS
# 0, 4, 5
REPRO_BRYO_MONO_P0 = """
    // find the megaspores
    if (individual.tag == 3) { //megaspores
        eggs = GAM_ARCHEGONIA_PER;
        //fertilize each egg
        for (rep in 1:eggs) {
            sperm = p0.sampleIndividuals(1, tag=2); //find a microspore
            //NOTE: each microspore gives rise to a male gametophyte, which will
            //produce many antheridia, giving rise to thousands of clonal sperm
            //because of this, sperm is not removed from the mating pool when used
            
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }
    else if (individual.tag == 1) { //already formed eggs
        sperm = p0.sampleIndividuals(1, tag=2); //find a microspore 
        if (sperm.size() == 1) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
            child.tag=0;
            if (GAM_MATERNAL_EFFECT > 0)
                child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }
    }
    else if (individual.tag == 4) { //add new gametophyte clones to p1 as haploids
        for (i in 1:GAM_CLONES_PER)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
        
        //also fertilize each egg
        eggs = GAM_ARCHEGONIA_PER;
        for (rep in 1:eggs) {
            sperm = p0.sampleIndividuals(1, tag=2); //find a microspore
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }
    else if (individual.tag == 6){ //performs gametophytic selfing
        eggs = GAM_ARCHEGONIA_PER;
        for (rep in 1:eggs){
        child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL,  NULL);
        child.tag=0;
        // Mother's fitness affects sporophyte fitness; see survival()
        if (GAM_MATERNAL_EFFECT > 0)
            child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }
    }

"""

LATE_BRYO_MONO = """
        //tag gametophytes that will clone
        gametophytes = c(p0.individuals[p0.individuals.tag==3], p0.individuals[p0.individuals.tag==2]);
        choose = rbinom(1, length(gametophytes), GAM_CLONE_RATE);
        clones = sample(gametophytes, choose);
        clones.tag = 4; //tag clones;
        
        //tag gametophytes that will self
        p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*GAM_SELF_RATE));
        clones.tag = 6; //tag selfed;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        //calculate chance a sporophyte will make a selfing egg
        choose = rbinom(1, p1.individualCount, SPO_SELF_CHANCE);
        selfed = p1.sampleIndividuals(choose);
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""
