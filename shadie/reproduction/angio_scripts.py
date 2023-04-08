#!/usr/bin/env python

"""
Spermatophyte (just angiosperm currently) string substitutions.
"""

#shadie DEFINITIONS
DEFS_ANGIO="""
// p0 = haploid population
// p1 = diploid population
// >1M- tmp parental tags
// <2000000 = female gametophyte (1N)
// >2000000 = male gametophyte (1N)
// 2 = gametophyte clones (1N) tag -------------- None in this model
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag
// 0 = used pollen (1N) tag
"""
#
# spo_pop_size
# spo_female_to_male_ratio
#
EARLY1_ANGIO_DIO = """
    sim.addSubpop('p1', SPO_POP_SIZE);
    sim.addSubpop('p0', 0);
    
    //set male and female flowers
    fem_num = asInteger(SPO_FEMALE_TO_MALE_RATIO*SPO_POP_SIZE);
    male_num = SPO_POP_SIZE-fem_num;
    
    spo_sex_starts = c(rep(1000001, fem_num), rep(2000001, male_num));
    for (idx in SPO_POP_SIZE)
        p1.individuals.tag = spo_sex_starts;
"""

EARLY_P0_FITNESS = """
    males = p0.individuals[p0.individuals.tag ==2];
    p0.fitnessScaling = (GAM_POP_SIZE / (p0.individualCount-length(males)));
    control = sample(males, asInteger(length(males)*POLLEN_CONTROL));
    control.fitnessScaling = 0.0;
"""

ANGIO_DIO_FITNESS_SCALE = """males = length(p0.individuals[p0.individuals.tag ==2]);
        p0.fitnessScaling = (gam_pop_size / (p0.individualCount-males));"""


#PARAMETERS
# -------------------------
# TAGS
# 2
ANGIO_P0_SURV = """
    if (individual.tag == 2)
        return T;
"""

# PARAMETERS
# spo_flowers_per
# flower_ovules_per
# flower_anthers_per
# anther_pollen_per
# spo_clones_per
#spo_maternal_effect
# -------------------------
# TAGS
# 1, 2, 41, 42
REPRO_ANGIO_DIO_P1 = """
    ind = individual;
    
    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < SPO_CLONE_RATE) {
        for (i in 1:SPO_CLONES_PER) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome2, NULL, NULL);
            child.tag = 4;
            child.setValue("parentid", ind.tag);
        }
    }
    
    
    if (individual.tag == 1000001)
    {
        // determine how many ovules were fertilized, out of the total
        meiosis_reps = rbinom(1, SPO_ARCHEGONIA_PER, FERTILIZATION_RATE);
        
        //one egg per archegonia (fertilized ovule)
        for (rep in meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            egg = p0.addRecombinant(ind.genome1, ind.genome2, breaks, NULL, NULL, NULL);
            egg.tag = 1000001;
        }
    
    }
    else //individual is male
    {
        meiosis_reps = floor(SPO_POLLEN_PER/4);
        for (rep in meiosis_reps)
        {
            breaks1 = sim.chromosome.drawBreakpoints(ind);
            breaks2 = sim.chromosome.drawBreakpoints(ind);
            
            // create four meiotic products
            child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL);
            child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
            child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
            child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
            children = c(child1, child2, child3, child4);
            children.tag = 2000001;
        }
    }
    if (SPO_MATERNAL_EFFECT > 0)
        children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
"""

# PARAMETERS
# pollen_per_stigma
#spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2
REPRO_ANGIO_DIO_P0  = """
    // iterate over each egg to find a mate (self, sib, or outcross)
    
    //Reproduction scripts run only females
    if (individual.tag < 2000000) {
    
    // get all males that could fertilize an egg of this female
    males = p0.individuals[p0.individuals.tag > 2000000];
        
        if (POLLEN_COMPETITION == T) {
            
            // sperm land on stigma
            pollen_pool = sample(males, POLLEN_PER_STIGMA);
            for (pollen in pollen_pool) {
                // store fitness value
                pollen.setValue("fitness", p0.cachedFitness(pollen.index));
                //pollen.tag = 0; do not remove for now
            }
            
            if (length(pollen_pool)>0) {
                //sort pollens by fitness
                fitness_vector = pollen_pool.getValue("fitness");
                sorted_fitness_vector = sort(fitness_vector, ascending=F);
                
                //calculate how many pollens attempt to fertilize
                attempts = 0;
                
                for (i in range(1:length(pollen_pool))) {
                    attempts = attempts + 1;
                    if (runif(1)<POLLEN_SUCCESS_RATE)
                        break;
                }
                idx = attempts-1;
                target_fitness = sorted_fitness_vector[idx];
                winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
                sperm = winners[0];
            }
        }
        
        else {//no pollen competition
            
            // get all males that could fertilize an egg of this female
            males = p0.individuals[p0.individuals.tag == 2000001];
            
            //only out-crossing is possible     
            // try at most 10 times to find a non-sib sperm, then skip.
            for (trial in 1:10) {
                        sperm = sample(males, 1);
                        child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                        child.tag = ifelse(runif(1)<FEMALE_TO_MALE_RATIO, 1000001, 2000001);
                    }
        }
    }
"""


# PARAMETERS
# spo_flowers_per
# flower_ovules_per
# flower_anthers_per
# anther_pollen_per
# spo_clones_per
# spo_maternal_effect
# -------------------------
# TAGS
# 0, 1, 2, 44, 5, 45
REPRO_ANGIO_MONO_P1="""
    ind = individual;
    
    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < SPO_CLONE_RATE) {
        for (i in 1:SPO_CLONES_PER) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome2, NULL, NULL);
            child.tag = 4;
            child.setValue("parentid", ind.tag);
        }
    }
    
    // parent tag is 1M + the parents unique index
    //female tags start with 1
    ftag = 1000000 + ind.index;
    //male tags start with 2
    mtag = 2000000 + ind.index;
    
    // determine how many ovules were fertilized, out of the total
    meiosis_reps = rbinom(1, SPO_ARCHEGONIA_PER, FERTILIZATION_RATE);
    
    //one egg per archegonia (fertilized ovule)
    for (rep in meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        egg = p0.addRecombinant(ind.genome1, ind.genome2, breaks, NULL, NULL, NULL);
        egg.tag = ftag;
    }
    
    
    meiosis_reps = floor(SPO_POLLEN_PER/4);
    for (rep in meiosis_reps)
    {
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        
        // create four meiotic products
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
        children = c(child1, child2, child3, child4);
        children.tag = mtag;
    }
    if (SPO_MATERNAL_EFFECT > 0)
        children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
"""

# PARAMETERS
# pollen_per_stigma
#spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2, 44, 5
REPRO_ANGIO_MONO_P0="""
    // iterate over each egg to find a mate (self, sib, or outcross)
    
    //Reproduction scripts run only females
    if (individual.tag < 2000000) {
    
    // get all males that could fertilize an egg of this female
    males = p0.individuals[p0.individuals.tag > 2000000];
        
        if (POLLEN_COMPETITION == T) {
            
            // sperm land on stigma
            pollen_pool = sample(males, POLLEN_PER_STIGMA);
            for (pollen in pollen_pool) {
                // store fitness value
                pollen.setValue("fitness", p0.cachedFitness(pollen.index));
                //pollen.tag = 0; do not remove for now
            }
            
            if (length(pollen_pool)>0) {
                //sort pollens by fitness
                fitness_vector = pollen_pool.getValue("fitness");
                sorted_fitness_vector = sort(fitness_vector, ascending=F);
                
                //calculate how many pollens attempt to fertilize
                attempts = 0;
                
                for (i in range(1:length(pollen_pool))) {
                    attempts = attempts + 1;
                    if (runif(1)<POLLEN_SUCCESS_RATE)
                        break;
                }
                idx = attempts-1;
                target_fitness = sorted_fitness_vector[idx];
                winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
                sperm = winners[0];
            }
        }
        
        else {//no pollen competition
            
            //only out-crossing is possible     
            // try at most 10 times to find a non-sib sperm, then skip.
            for (trial in 1:10) {
                sperm = sample(males, 1);
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag = 3;
            }
        }
    }
"""