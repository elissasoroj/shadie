#!/usr/bin/env python

"""
Spermatophyte (just angiosperm currently) string substitutions.
"""

#
# spo_pop_size
# spo_female_to_male_ratio
#
EARLY1_ANGIO = """
    // diploid sporophyte pop
    sim.addSubpop('p1', spo_pop_size);

    // haploid gametophyte pop
    sim.addSubpop('p0', 0);

    // tag individuals as male or female.
    fems = spo_female_to_male_ratio * spo_pop_size;
    spo_sex_starts = c(rep(1, asInteger(fems)), 
        rep(0, asInteger(spo_pop_size-fems)));
    p1.individuals.tag = spo_sex_starts;
"""

# PARAMETERS
# spo_flowers_per
# flower_ovules_per
# flower_anthers_per
# anther_microspores_per
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
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_microspores_per/4);
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
        
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_microspores_per/4);
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
# anther_microspores_per
# spo_clones_per
# spo_maternal_effect
# -------------------------
# TAGS
# 0, 1, 2, 44, 5, 45
REPRO_ANGIO_MONO_P1="""
    g_1 = genome1;
    g_2 = genome2;
    //sporophytes make hermaphroditic flowers
    if (individual.tag == 0) {
        //make ovules
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

        //make microspores (pollen)
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_microspores_per/4);
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
    
    if (individual.tag == 44) { //save cloned spo to p0
        //make the clones
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 44;
        
        //clone can reproduce normally as well
        //make ovules
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

        //make microspores (pollen)
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_microspores_per/4);
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

    //make sporophytic selfed
    if (individual.tag == 5) {
        // for each ovule generated by the sporophyte, perform meiosis twice
        meiosis_reps = asInteger(spo_flowers_per*flower_ovules_per);
        for (rep in 1:meiosis_reps) {

            // sample meiosis crossover position to generate 4 microspores
            // male outcross
            breaks_m = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks_m).tag = 2;
            child2 = p0.addRecombinant(NULL, NULL, NULL, genome_1, genome_2, breaks_m).tag = 2;
            child3 = p0.addRecombinant(genome_2, genome_1, breaks_m, NULL, NULL, NULL).tag = 2;

            //only one megaspore will be produced, used for the new selfed sporophyte
            breaks_f = sim.chromosome.drawBreakpoints(individual);
            
            // add the diploid selfed 
            p0.addRecombinant(genome_1, genome_2, breaks_m, genome_2, genome_1, breaks_f).tag = 5;

            children = c(child1, child2, child3);
            children.tag = 2;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }

        #make remaining microspores
        //make microspores (pollen)
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_microspores_per/4);
        male_meiosis_reps = microsporocytes-meiosis_reps
        //perform meiosis for each microsporocyte to produce microspores, which will mature into pollen
        for (rep in 1:male_meiosis_reps){
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

    if (individual.tag == 45) {
         //make the clones
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 44;

        //make selfed
        // for each ovule generated by the sporophyte, perform meiosis twice
        meiosis_reps = asInteger(spo_flowers_per*flower_ovules_per);
        for (rep in 1:meiosis_reps) {

            // sample meiosis crossover position to generate 4 microspores
            // male outcross
            breaks_m = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks_m).tag = 2;
            child2 = p0.addRecombinant(NULL, NULL, NULL, genome_1, genome_2, breaks_m).tag = 2;
            child3 = p0.addRecombinant(genome_2, genome_1, breaks_m, NULL, NULL, NULL).tag = 2;

            //only one megaspore will be produced, used for the new selfed sporophyte
            breaks_f = sim.chromosome.drawBreakpoints(individual);
            
            // add the diploid selfed 
            p0.addRecombinant(genome_1, genome_2, breaks_m, genome_2, genome_1, breaks_f).tag = 5;

            children = c(child1, child2, child3);
            children.tag = 2;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }

        #make remaining microspores
        //make microspores (pollen)
        microsporocytes = asInteger(spo_flowers_per*flower_anthers_per*anther_microspores_per/4);
        male_meiosis_reps = microsporocytes-meiosis_reps
        //perform meiosis for each microsporocyte to produce microspores, which will mature into pollen
        for (rep in 1:male_meiosis_reps){
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
# 1, 2, 44, 5
REPRO_ANGIO_MONO_P0="""
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
            p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL).tag=0;
            sperm.tag = 20; //sperm goes into used pool
        }
    }

    //move clones back into p1,reset tag
    if (individual.tag == 44)
        p1.addCloned(individual).tag=0

    //move clones back into p1, reset tag
    if (individual.tag == 5)
        p1.addCloned(individual).tag=0

"""

# PARAMETERS
# spo_clone_rate
# spo_self_rate
# -------------------------
# TAGS
# 44, 5, 45
LATE_ANGIO_MONO = """
}
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        clones.tag == 44;
        
        //sporophytic selfing
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds
    }
"""

