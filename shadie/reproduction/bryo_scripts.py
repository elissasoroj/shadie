#!/usr/bin/env python

"""
Bryophyte specific SLIM script snippets used for string substitution.
"""
S4_TAG = """individual.tag = ifelse(runif(1)<0.5, 3, 2);"""

FUNCTIONS_BRYO_DIO = """
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
        clones.tag = 4; //tag clones;
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
    // normal sporophyte: make both female and male spores
    if (individual.tag == 0) { //focus on one sporophyte

        // each meiotic division produces 4 spores
        meiosis_reps = asInteger(spo_spores_per/4);
        for (rep in 1:meiosis_reps) {

            // sample a meiosis crossover position.
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create 2 megaspores and 2 microspores (all recombinant for now)
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 0;
        }
    }

    // cloned gametophytes from last gen: copy and move clones directly to p0,
    //these cloned gametophyte becomes a 'normal' gametophyte in p0, tag=4 removed
    if (individual.tag == 4) 
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 0;

    // selfing sporophyte: make selfed each meiotic division
    if (individual.tag == 5) {
        meiosis_reps = asInteger(spo_spores_per/4);
        for (rep in 1:meiosis_reps)
        breaks1 = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome1, genome2, breaks1, NULL, NULL, NULL).tag = 0;p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 0;

        breaks2 = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(genome1, genome2, breaks2, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome1, genome2, breaks2, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome2, genome1, breaks2, NULL, NULL, NULL).tag = 0;

        // add the diploid selfed 
        p0.addRecombinant(genome1, genome2, breaks1, genome2, genome1, breaks2).tag = 5;

        }
"""

## PARAMETERS
# 
# ------------------
# TAGS
# 0, 4, 5
REPRO_BRYO_MONO_P0 = """
    // p0 at this point contains microspores and megaspores
    // normal gametophyte: makes both archegonia and antheridia
    if (individual.tag == 0) {

        // each megaspore makes N archegonia, each of which contains one egg
        eggs = gam_archegonia_per;
        for (egg in 1:eggs) {

            // sample a sperm from p0; each microspore creates many sperm
            sperm = p0.sampleIndividuals(1);
            //note, the sperm will NOT be removed from the pool, because it is actually
            //a pool of microspores, which can produce thousands of clonal sperm
            
            if (sperm.size() == 1) {

                // create zygote and set tag to 0
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag = 0;

                // store maternal fitness to apply maternal effects later (in survival).
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
    }

    // gametophytic clone: add gametophyte to p1 and set tag=4.
    if (individual.tag == 4) {
        for (i in 1:gam_clones_per) {
            p1.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
        }
    }

    // sporophytic selfed: move into p1 and set tag=0.
    if (individual.tag == 5) 
        p1.addCloned(individual).tag = 0;
"""

LATE_BRYO_MONO = """
//odd = starts with gam in p0, generates spo into p1
    p0_size = length(p0.individuals);
    //tag gametophytic clones
    p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate)).tag=4;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        number_selfed = rbinom(1, length(p1_size), spo_self_chance);
        selfed = p1.sampleIndividuals(number_selfed);
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""
