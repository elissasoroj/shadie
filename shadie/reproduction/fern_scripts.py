#!/usr/bin/env python

"""
Pteridophyte specific SLIM script snippets used for string substitution.
"""
PTER_HETERO_FITNESS_SCALE = """males = length(p0.individuals[p0.individuals.tag ==2]);
        p0.fitnessScaling = (GAM_POP_SIZE / (p0.individualCount-males));"""

# PARAMETERS
#spo_spore_per
#spo_maternal_effect
#gam_archegonia_per
# -------------------------
# TAGS
# 0, 5
FUNCTIONS_PTER_HOMOSPORE = """
    //p0 = haploid population
    //p1 = diploid population

    //0 = hermaphrodite
    //1 = female
    //2 = male
    //4 = gametophyte clone
    //44 = sporophyte clone
    //45 = sporophyte cloned and selfed
    //5 = sporophytic selfed
    //6 = gametophyte selfed
   
   //shadie-ddefined functions
   // generates gametes from sporophytes
function (void)generate_gametes(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
        
        children = c(child1, child2, child3, child4);
        children.tag = 0;
        
        // Mother's fitness affects gametophyte fitness; see survival()
        if (SPO_MATERNAL_EFFECT > 0)
            children.setValue("maternal_fitness", ind.subpopulation.cachedFitness(ind.index));
    }
}

function (void)sporophyte_selfs(object<Individual>$ ind){
    breaks_m1 = sim.chromosome.drawBreakpoints(ind);
    breaks_m2 = sim.chromosome.drawBreakpoints(ind);
    breaks_f1 = sim.chromosome.drawBreakpoints(ind);
    breaks_f2 = sim.chromosome.drawBreakpoints(ind);
    //sporophyte makes a spore
    for (i in 1:SPO_SPORES_PER){
        //spore is made - germinates into gametophyte
        for (i in 1:GAM_ARCHEGONIA_PER){ //gametophyte makes archegonia
            // 4 gametophytes that produce archegonia+antheridia per meiosis
            child1 = p0.addRecombinant(ind.genome2, ind.genome1, breaks_m1, NULL, NULL, NULL);
            child2 = p0.addRecombinant(ind.genome1, ind.genome2, breaks_m1, NULL, NULL, NULL);
            child3 = p0.addRecombinant(ind.genome2, ind.genome1, breaks_m2, NULL, NULL, NULL);
            
            //4 gametophytes that produce archegonia+antheridia per meiosis
            child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks_f1, NULL, NULL, NULL);
            child5 = p0.addRecombinant(ind.genome1, ind.genome2, breaks_f1, NULL, NULL, NULL);
            child6 = p0.addRecombinant(ind.genome2, ind.genome1, breaks_f2, NULL, NULL, NULL);
            
            p1.addRecombinant(ind.genome1, ind.genome2, breaks_m2, ind.genome1, ind.genome2, breaks_f2).tag = 5;
            
            //save all children produced to vector
            children = c(child1, child2, child3, child4, child5, child6);
            children.tag = 0;
            
            // Mother's fitness affects gametophyte fitness; see survival()
            if (SPO_MATERNAL_EFFECT > 0)
                children.setValue("maternal_fitness", ind.subpopulation.cachedFitness(individual.index));
        }
    }
}
"""

# PARAMETERS
#spo_spore_per
#spo_maternal_effect
#spo_clones_per
#gam_archegonia_per
# -------------------------
# TAGS
# 0, 4, 44, 5, 45, 6, 46
REPRO_PTER_HOMOSPORE_P1 = """
    if (individual.tag == 0)
        generate_gametes(individual, asInteger(SPO_SPORES_PER/4));
    
    else if (individual.tag == 44) { //save cloned spo to p0
        //make the clones
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        
        //individual can make gametes too
        generate_gametes(individual, asInteger(SPO_SPORES_PER/4));
    }
    
    //sporophytic selfing
    else if (individual.tag == 5) {
        //calculate how many spores will 
        //sporophyte selfs
        sporophyte_selfs(individual);
        
        //individual undergoes remaining meiosis
        generate_gametes(individual, asInteger(SPO_SPORES_PER/4) - 2);
    }
    
    else if (individual.tag == 45) {
        //make clones and move to p0 as diploids
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        
        //individual selfs
        sporophyte_selfs(individual);
        
        //individual undergoes remaining meiosis
        generate_gametes(individual, asInteger(SPO_SPORES_PER/4) - 2);
    }
"""

# PARAMETERS
#gam_archegonia_per
#gam_clones_per
# -------------------------
# TAGS
# 0, 1, 2, 4, 5, 6
REPRO_PTER_HOMOSPORE_P0 = """
    //all gametophytes are hermaphroditic
    if (individual.tag == 0) {
        //for each egg per gametophyte, perform fertilization
        for (rep in 1:GAM_ARCHEGONIA_PER){
            sperm = p0.sampleIndividuals(1); //all individuals make sperm;
            //NOTE: each gametophyte makes many antheridia, each of which 
            //produec thousands of clonal sperm. For this reason, sperm is 
            //not removed from the mating pool once used
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                
                // Mother's fitness affects sporophyte fitness; see survival()
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", individual.subpopulation.cachedFitness(individual.index));
            }
        }
    }
    
    else if (individual.tag == 4) { //add new gametophyte clones to p1 as haploids
        for (i in 1:GAM_CLONES_PER)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }
    
    else if (individual.tag == 6){ //performs gametophytic selfing
        child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL,  NULL);
        child.tag=0;
        // Mother's fitness affects sporophyte fitness; see survival()
        if (GAM_MATERNAL_EFFECT > 0)
            child.setValue("maternal_fitness", individual.subpopopulation.cachedFitness(individual.index));
    }
"""

# PARAMETERS
#gam_archegonia_per
#spo_clone_rate
#gam_clone_rate
#spo_self_rate
#gam_self_rate
# -------------------------
# TAGS
# 0, 1, 2, 4, 44, 5, 45, 6, 46
LATE_PTER_HOMOSPORE = """
    //tag gametophytes that will clone
        p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*GAM_CLONE_RATE));
        clones.tag = 4; //tag clones;
        
        //tag gametophytes that will self
        p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*GAM_SELF_RATE));
        clones.tag = 6; //tag selfed;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        //tag sporophytes that will clone
        clones = p1.sampleIndividuals(asInteger(p1_size*SPO_CLONE_RATE));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo;
        
        //tag sporophytes that will self
        number_selfed = rbinom(1, length(p1_size), SPO_SELF_RATE);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds;
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""

#-----------------------------------------------------------------------
FUNCTIONS_PTER_HETEROSPORE = """
//p0 = haploid population
//p1 = diploid population

//0 = hermaphrodite
//1 = female
//2 = male
//4 = gametophyte clone
//44 = sporophyte clone
//45 = sporophyte cloned and selfed
//5 = sporophytic selfed
//6 = gametophyte selfed

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
    strobilus_female_ratio = RS_MEGASPORANGIA_PER/(RS_MICROSPORANGIA_PER+RS_MEGASPORANGIA_PER);
    count = 0;
    if (runif(1) < strobilus_female_ratio){
        for (i in 1:MEGASPORANGIA_MEGASPORES_PER){
            for (i in 1:GAM_ARCHEGONIA_PER){ //gametophyte makes archegonia
                if (runif(1)<SPO_SELF_RATE){
                    count = count + 1;
                    breaks1 = sim.chromosome.drawBreakpoints(individual);
                    breaks2 = sim.chromosome.drawBreakpoints(individual);
                    breaks_f = sim.chromosome.drawBreakpoints(individual);
                    
                    //4 gametophytes that produce antheridia per meiosis
                    child1 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
                    child2 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
                    child3 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);
                    
                    p1.addRecombinant(ind.genome1, ind.genome2, breaks1, ind.genome1, ind.genome2, breaks_f).tag = 5;
                    
                    //save all children produced to vector
                    children = c(child1, child2, child3);
                    children.tag = 0;
                    
                    // Mother's fitness affects gametophyte fitness; see survival()
                    if (SPO_MATERNAL_EFFECT > 0)
                        children.setValue("maternal_fitness", ind.subpopulation.cachedFitness(individual.index));
                }
                else
                    make_eggs(individual, 1);
            }
        }
    }
    else{
        meiosis_reps = asInteger(MICROSPORANGIA_MICROSPORES_PER*RS_MICROSPORANGIA_PER)-count;
        make_microspores(individual, meiosis_reps);
    }
}
"""


# PARAMETERS
#rs_megasporangia_per
#rs_microsporangia_perv
#megasporangia_megaspores_per
#microsporangia_microspores_per
#spo_clones_per
# -------------------------
# TAGS
# 0, 1, 2, 44, 5, 45
REPRO_PTER_HETEROSPORE_P1 = """
    if (individual.tag == 0){
        strobilus_female_ratio = RS_MEGASPORANGIA_PER/(RS_MICROSPORANGIA_PER+RS_MEGASPORANGIA_PER);
        if (runif(1) < strobilus_female_ratio){
            meiosis_reps = asInteger(MEGASPORANGIA_MEGASPORES_PER*RS_MEGASPORANGIA_PER);
            make_eggs(individual, meiosis_reps);
        }
        else{
            meiosis_reps = asInteger((RS_MICROSPORANGIA_PER*MICROSPORANGIA_MICROSPORES_PER)/4);
            make_microspores(individual, meiosis_reps);
        }
    }
    
    if (individual.tag == 44) { //save cloned spo to p0
        //sporophyte clones
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        
        //cloned individual can make gametes too
        strobilus_female_ratio = RS_MEGASPORANGIA_PER/(RS_MICROSPORANGIA_PER+RS_MEGASPORANGIA_PER);
        if (runif(1) < strobilus_female_ratio){
            meiosis_reps = asInteger(MEGASPORANGIA_MEGASPORES_PER*RS_MEGASPORANGIA_PER);
            make_eggs(individual, meiosis_reps);
        }
        else{
            meiosis_reps = asInteger((RS_MICROSPORANGIA_PER*MICROSPORANGIA_MICROSPORES_PER)/4);
            make_microspores(individual, meiosis_reps);
        }
    }
    
    //sporophytic selfing
    if (individual.tag == 5) //sporophyte selfs
        sporophyte_selfs(individual);
    
    if (individual.tag == 45) {
        //make clones and move to p0 as diploids
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        
        //sporophyte selfs
        sporophyte_selfs(individual);
    
    }
"""

# PARAMETERS
#gam_archegonia per (female gam = mature megaspore)
# -------------------------
# TAGS
#0, 1, 2, 44, 5, 20
REPRO_PTER_HETEROSPORE_P0 = """
    // find the megaspores
    if (individual.tag == 1) {
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
            
            }
        }
    }
    else if (individual.tag == 4) { //add new gametophyte clones to p1 as haploids
        for (i in 1:GAM_CLONES_PER)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }
"""

# PARAMETERS
#spo_clone_rate
#spo_self_rate
# -------------------------
# TAGS
# 44, 45, 5
LATE_PTER_HETEROSPORE = """
    //tag gametophytes that will clone
        p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*GAM_CLONE_RATE));
        clones.tag = 4; //tag clones;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        //tag sporophytes that will clone
        clones = p1.sampleIndividuals(asInteger(p1_size*SPO_CLONE_RATE));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo;
        
        //tag sporophytes that will self
        number_selfed = rbinom(1, length(p1_size), SPO_SELF_RATE);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds;
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""