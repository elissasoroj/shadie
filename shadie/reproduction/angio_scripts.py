#!/usr/bin/env python

"""
Spermatophyte (just angiosperm currently) string substitutions.
"""

FUNCTIONS_ANGIO_MONO="""
//p0 = haploid population
//p1 = diploid population

//0 = hermaphrodite
//1 = female
//2 = male
//3 = megaspore
//44 = sporophyte clone
//45 = sporophyte cloned and selfed
//5 = sporophytic selfed

//shadie-defined functions
// generates gametes from sporophytes
function (void)make_microspores(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){
        //4 microspores per meiosis rep
        breaks1 = sim.chromosome.drawBreakpoints(individual);
        breaks2 = sim.chromosome.drawBreakpoints(individual);
        child1 = p0.addRecombinant(individual.genome1, individual.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(individual.genome2, individual.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(individual.genome1, individual.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(individual.genome2, individual.genome1, breaks2, NULL, NULL, NULL);
        
        children = c(child1, child2, child3, child4);
        children.tag = 2;
        
        // Mother's fitness affects gametophyte fitness; see survival()
        if (SPO_MATERNAL_EFFECT > 0)
            children.setValue("maternal_fitness", individual.subpopulation.cachedFitness(individual.index));
    
    }
}

function (void)make_eggs(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){
        breaks = sim.chromosome.drawBreakpoints(individual);
        child1 = p0.addRecombinant(individual.genome1, individual.genome2, breaks, NULL, NULL, NULL);
        child1.tag = 1;
        
        // Mother's fitness affects gametophyte fitness; see survival()
        if (SPO_MATERNAL_EFFECT > 0)
            child1.setValue("maternal_fitness", individual.subpopulation.cachedFitness(individual.index));
    }
}

function (void)sporophyte_selfs(object<Individual>$ ind){
    //sporophyte makes a megasporangia
    eggs = FLOWER_OVULES_PER*SPO_FLOWERS_PER;
    eggs_selfed = eggs*EGG_SPO_SELF_RATE;
    for (i in 1:eggs_selfed){
        breaks1 = sim.chromosome.drawBreakpoints(individual);
        breaks2 = sim.chromosome.drawBreakpoints(individual);
        breaks_f = sim.chromosome.drawBreakpoints(individual);
        
        //4 microspores produced
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
        
        //only one egg produced, which selfs
        p1.addRecombinant(ind.genome1, ind.genome2, breaks1, ind.genome1, ind.genome2, breaks_f).tag = 5;
        
        //save all children produced to vector
        children = c(child1, child2, child3, child4);
        children.tag = 2;
        
        // Mother's fitness affects gametophyte fitness; see survival()
        if (SPO_MATERNAL_EFFECT > 0)
            children.setValue("maternal_fitness", ind.subpopulation.cachedFitness(individual.index));
    }
    //make non-selfed eggs
    make_eggs(individual, eggs-eggs_selfed);
    
    //make the rest of the pollen
    meiosis_reps = asInteger(SPO_FLOWERS_PER*FLOWER_ANTHERS_PER*ANTHER_POLLEN_PER/4) - eggs_selfed;
    make_microspores(individual, meiosis_reps);
}
"""
#
# spo_pop_size
# spo_female_to_male_ratio
#
EARLY1_ANGIO = """
    sim.addSubpop('p1', spo_pop_size);
    sim.addSubpop('p0', 0);

    // tag individuals as male or female.
    fems = spo_female_to_male_ratio * spo_pop_size;
    spo_sex_starts = c(rep(1, asInteger(fems)), 
        rep(0, asInteger(spo_pop_size-fems)));
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
    g_1 = genome1;
    g_2 = genome2;
    //female sporophyte makes megaspores (one ovule = one megaspore = one egg)
    if (individual.tag == 1) {
        meiosis_reps = asInteger(spo_flowers_per*flower_ovules_per);
        //each ovule undergoes meiosis once and produces a single megaspore
        for (rep in 1:meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
            child1.tag = 1;
            
            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }
    
    //male sporophytes produce microspores (one microsporocyte -> 4 microspores = 4 pollen)
    if (individual.tag == 2) {
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_pollen_per/4);
        //perform meiosis for each microsporocyte to produce microspores, which will mature into pollen
        for (rep in 1:microsporocytes){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4);
            children.tag = 2;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }
    
    if (individual.tag == 41) { //save cloned spo to p0
        //make the clones
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 41;
        
        meiosis_reps = asInteger(spo_flowers_per*flower_ovules_per);
        //each ovule undergoes meiosis once and produces a single megaspore
        for (rep in 1:meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
            child1.tag = 1;
            
            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }
    
    //else made a male-only gametophyte
    if (individual.tag == 42) {
        //make the clones
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 42;
        
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_pollen_per/4);
        //perform meiosis for each microsporocyte to produce microspores, which will mature into pollen
        for (rep in 1:microsporocytes){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4);
            children.tag = 2;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }
"""

# PARAMETERS
# pollen_per_stigma
#spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2
REPRO_ANGIO_DIO_P0  = """
    // females find male gametes to reproduce
    if (individual.tag == 1) {
        if (pollen_comp == T) {

            // sperm land on stigma
            pollen_pool = p0.sampleIndividuals(pollen_per_stigma, tag=2);
            for (pollen in pollen_pool) {
                // store fitness value
                pollen.setValue("fitness", p0.cachedFitness(pollen.index));
                pollen.tag = 2;
            }

            if (length(pollen_pool)>0) {
                //sort pollens by fitness
                fitness_vector = pollen_pool.getValue("fitness");
                sorted_fitness_vector = sort(fitness_vector, ascending=F);
                                
                //calculate how many pollens attempt to fertilize
                attempts = 0;

                for (i in range(1:length(pollen_pool))) {
                    attempts = attempts + 1;
                    if (runif(1)<pollen_success_rate)
                        break;
                }
                idx = attempts-1;    
                target_fitness = sorted_fitness_vector[idx];
                winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
                sperm = winners[0];
            }
        }

        else //no pollen competition
            sperm = p0.sampleIndividuals(1, tag=2); //find pollen to fetilize ovule

        if (sperm.size() == 1) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
            sperm.tag = 20; //sperm goes into used pool

            if (runif(1) <= spo_female_to_male_ratio)
                child.tag = 1;
            else
                child.tag = 0;
        }
    }
"""

# PARAMETERS
# spo_clone_rate
# -------------------------
# TAGS
# 41, 42
LATE_ANGIO_DIO = """
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        femclones = clones[clones.tag == 1];
        femclones.tag = 41;
        maleclones = clones[clones.tag == 2];
        maleclones.tag = 42;
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
    if (individual.tag == 0){
        f_meiosis_reps = FLOWER_OVULES_PER*SPO_FLOWERS_PER;
        make_eggs(individual, f_meiosis_reps);
        
        m_meiosis_reps = asInteger(SPO_FLOWERS_PER*FLOWER_ANTHERS_PER*ANTHER_POLLEN_PER/4);
        make_microspores(individual, m_meiosis_reps);
    }
    
    else if (individual.tag == 44) { //save cloned spo to p1
        //sporophyte clones
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        
        //cloned individual can make gametes too
        f_meiosis_reps = FLOWER_OVULES_PER*SPO_FLOWERS_PER;
        make_eggs(individual, f_meiosis_reps);
        
        m_meiosis_reps = asInteger(SPO_FLOWERS_PER*FLOWER_ANTHERS_PER*ANTHER_POLLEN_PER/4);
        make_microspores(individual, m_meiosis_reps);
    }
    
    //sporophytic selfing
    else if (individual.tag == 5) //sporophyte selfs
        sporophyte_selfs(individual);
    
    else if (individual.tag == 45) {
        //make clones and move to p0 as diploids
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        
        //sporophyte selfs
        sporophyte_selfs(individual);
    }
"""

# PARAMETERS
# pollen_per_stigma
#spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2, 44, 5
REPRO_ANGIO_MONO_P0="""
    // find an egg
    if (individual.tag == 1) {
        if (POLLEN_COMP == T) {
            
            // sperm land on stigma
            pollen_pool = p0.sampleIndividuals(pollen_comp_stigma_pollen_per, tag=2);
            for (pollen in pollen_pool) {
                // store fitness value
                pollen.setValue("fitness", p0.cachedFitness(pollen.index));
                pollen.tag = 20;
            }
            
            if (length(pollen_pool)>0) {
                //sort pollens by fitness
                fitness_vector = pollen_pool.getValue("fitness");
                sorted_fitness_vector = sort(fitness_vector, ascending=F);
                
                //calculate how many pollens attempt to fertilize
                attempts = 0;
                
                for (i in range(1:length(pollen_pool))) {
                    attempts = attempts + 1;
                    if (runif(1)<pollen_success_rate)
                        break;
                }
                idx = attempts-1;
                target_fitness = sorted_fitness_vector[idx];
                winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
                sperm = winners[0];
            }
        }
        else{
            sperm = p0.sampleIndividuals(1, tag=2); //find a pollen
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                sperm.tag=20;
            }
        }
    }
"""

# PARAMETERS
# spo_clone_rate
# spo_self_chance
# -------------------------
# TAGS
# 44, 5, 45
LATE_ANGIO_MONO = """
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        //tag sporophytes that will clone
        clones = p1.sampleIndividuals(asInteger(p1_size*SPO_CLONE_RATE));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo;
        
        //tag sporophytes that will self
        number_selfed = rbinom(1, length(p1_size), SPO_SELF_CHANCE);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds;
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""

