#!/usr/bin/env python

"""
Bryophyte specific SLIM script snippets used for string substitution.
"""
S4_TAG = """individual.tag = ifelse(runif(1)<0.5, 3, 2);"""

FUNCTIONS_BRYO = """
// shadie DEFINITIONS
//p0 = haploid population
//p1 = diploid population

//0 = hermaphrodite
//1 = egg
//2 = micropore (male gametophyte)
//3 = megaspore (female gametophyte)
//4 = gametophyte clone
//5 = sporophytic selfed


function (void)make_spores(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){
        //4 microspores per meiosis rep
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL).tag=3; //megaspore
        p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL).tag=2;
        p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL).tag=3; //megaspore
        p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL).tag=2;
    }
}

function (void)sporophyte_selfs(object<Individual>$ ind){
    meiosis_reps = SPO_SPORES_PER/4;
    additional_reps = 0;
    for (rep in 1:meiosis_reps);{//for each meiosis rep
        if (runif(1) < 4*GAM_ARCHEGONIA_PER*EGG_SPO_SELF_RATE){//chance megaspoore will self
            
            eggs_made = 4*GAM_ARCHEGONIA_PER; //number of eggs produced by sporophyte that selfed
            //eggs_selfed = rbinom(1, eggs_made, EGG_SPO_SELF_RATE); //number of eggs selfed 
            megaspore1_eggs = rbinom(1, eggs_made, 1/4);
            megaspore2_eggs = rbinom(1, (eggs_made-megaspore1_eggs), 1/3);
            megaspore3_eggs = rbinom(1, eggs_made-(megaspore1_eggs+megaspore2_eggs), 1/2);
            megaspore4_eggs = eggs_made - (megaspore1_eggs+megaspore2_eggs+megaspore3_eggs);
            
            megaspore_eggs = c(megaspore1_eggs, megaspore2_eggs, megaspore3_eggs, megaspore4_eggs);
            
            selfed_eggs = c();
            for (i in 0:3)
                selfed_eggs = c(selfed_eggs, rbinom(1, megaspore_eggs[i], EGG_SPO_SELF_RATE)); //choose which eggs self
            
            breaks1 = sim.chromosome.drawBreakpoints(ind);
            breaks2 = sim.chromosome.drawBreakpoints(ind);
            breaks3 = sim.chromosome.drawBreakpoints(ind);
            breaks4 = sim.chromosome.drawBreakpoints(ind);
            
            //add all microspores to mating pool
            //p0.addRecombinant(ind.genome2, ind.genome1, breaks1).tag=1
            p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL).tag=2;
            //p0.addRecombinant(ind.genome2, ind.genome1, breaks2).tag=1
            p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL).tag=2;
            
            //p0.addRecombinant(ind.genome2, ind.genome1, breaks3).tag=1
            p0.addRecombinant(ind.genome1, ind.genome2, breaks3, NULL, NULL, NULL).tag=2;
            //p0.addRecombinant(ind.genome2, ind.genome1, breaks4).tag=1
            p0.addRecombinant(ind.genome1, ind.genome2, breaks4, NULL, NULL, NULL).tag=2;
            
            for (i in 0:3){ //for each megaspore
                for (egg in 1:(megaspore_eggs[i]-selfed_eggs[i])){ //make non-selfed eggs
                    if (i == 0)
                        breaks = breaks1;
                    else if (i == 1)
                        breaks = breaks2;
                    else if (i == 2)
                        breaks = breaks3;
                    else if (i == 3)
                        breaks = breaks4;
                    
                    p0.addRecombinant(ind.genome1, ind.genome2, breaks, NULL, NULL, NULL).tag=1;
                }
                for (egg in 1:selfed_eggs[i]){ //for each egg they make
                    if (i<2){
                        if (rbinom(1, 1, 0.5)>0)
                            fbreaks = breaks1;
                        else
                            fbreaks = breaks2;
                        
                        if (rbinom(1, 1, 0.5)>0)
                            mbreaks = breaks3;
                        else
                            mbreaks = breaks4;
                        p1.addRecombinant(ind.genome1, ind.genome2, fbreaks,
                            ind.genome2, ind.genome1, mbreaks).tag = 5;
                    }
                    else{
                        if (rbinom(1, 1, 0.5)>0)
                            fbreaks = breaks3;
                        else
                            fbreaks = breaks4;
                        
                        if (rbinom(1, 1, 0.5)>0)
                            mbreaks = breaks1;
                        else
                            mbreaks = breaks2;
                        p1.addRecombinant(ind.genome1, ind.genome2, fbreaks,
                            ind.genome2, ind.genome1, mbreaks).tag = 5;
                    }
                }
            }
        }
        else
            additional_reps = 1+additional_reps; // reps of meiosis that still need to be performed
    }
    make_spores(individual, additional_reps);
}

"""
# gam_female_to_male_ratio
# spo_spores_per
#gam_archegonia_per
#gam_antheridia_per
# -------------------------
# TAGS
# 0, 1, 2, 4, 5
REPRO_BRYO_DIO_P1 = """
    if (individual.tag == 0){
        meiosis_reps = asInteger(SPO_SPORES_PER/4);
        make_spores(individual, meiosis_reps);
    }
    
    //sporophytic selfing
    if (individual.tag == 5) //sporophyte selfs
        sporophyte_selfs(individual);
"""

# Parameters:
# gam_archegonia_per (we could use "megaspore_archegonia_per" or "gam_eggs_per")
# gam_clones_per
# -------------------------
# TAGS:
# 0, 1, 4, 5, 20
REPRO_BRYO_DIO_P0 = """
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
    else if (individual.tag == 41) { //add new gametophyte clones to p1 as haploids
        for (i in 1:GAM_CLONES_PER)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 41;
        
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
    else if (individual.tag == 42) { //add new gametophyte clones to p1 as haploids
        for (i in 1:GAM_CLONES_PER)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 42;
    }
"""

# PARAMETERS
# gam_clone_rate
# spo_self_chance
# ------------------
# TAGS
# 4, 5
LATE_BRYO_DIO = """
        //tag gametophytes that will clone
        gametophytes = c(p0.individuals[p0.individuals.tag==3], p0.individuals[p0.individuals.tag==2]);
        choose = rbinom(1, length(gametophytes), GAM_CLONE_RATE);
        clones = sample(gametophytes, choose);
        femclones = clones[clones.tag==3];
        maleclones=clones[clones.tag==2];
        femclones.tag = 41; //tag clones;
        maleclones.tag=42;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        //calculate chance a sporophyte will make a selfing egg
        p1_size = length(p1.individuals);
        choose = rbinom(1, p1_size, SPO_SELF_CHANCE);
        selfed = p1.sampleIndividuals(choose);
        selfed.tag = 5; //tag sporophytic selfing inds;
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
