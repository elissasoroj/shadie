#!/usr/bin/env python

"""
Pteridophyte specific SLIM script snippets used for string substitution.
"""

# PARAMETERS
#spo_spore_per
#spo_maternal_effect
#spo_clones_per
#gam_archegonia_per
# -------------------------
# TAGS
# 0, 1, 2, 4, 5
REPRO_PTER_HOMOSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    //normal sporophyte makes spores
    if (individual.tag == 0) {
    	meiosis_reps = asInteger(spo_spores_per/4)
    	for (rep in 1:reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4)
            children.tag = 0

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
    	}
    }
    
    if (individual.tag == 4) //move gam clones directly to p0
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 0;
    
    if (individual.tag == 44) { //save cloned spo to p0
        //make the clones
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 44;
        
        //individual can make gametes too
        if (individual.tag == 0) {
	    	meiosis_reps = asInteger(spo_spores_per/4)
	    	for (rep in 1:reps){
	            breaks = sim.chromosome.drawBreakpoints(individual);
	            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            
	            children = c(child1, child2, child3, child4)
	            children.tag = 0

	            // Mother's fitness affects gametophyte fitness; see survival()
	            if (spo_maternal_effect > 0){
	                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
	            }
	    	}
	    }
    }
    
    //sporophytic selfing
    // selfing sporophyte: make male, female, and self-recombinant spores
    if (individual.tag == 5) {

        //perform meiosis twice for the selfed diploid

        // 4 gametophytes that producearchegonia+antheridia per meiosis
        breaks_m = sim.chromosome.drawBreakpoints(individual);
        child5 = p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks_m)
        child6 = p0.addRecombinant(NULL, NULL, NULL, genome_1, genome_2, breaks_m)
        child7 = p0.addRecombinant(genome_2, genome_1, breaks_m, NULL, NULL, NULL)

        //4 gametophytes that produce archegonia+antheridia per meiosis
        breaks_f = sim.chromosome.drawBreakpoints(individual);
        child8 = p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks_m)
        child9 = p0.addRecombinant(NULL, NULL, NULL, genome_1, genome_2, breaks_m)
        child10 = p0.addRecombinant(genome_2, genome_1, breaks_m, NULL, NULL, NULL)

        // one gametophyte from each set above produces selfed
        p0.addRecombinant(genome_1, genome_2, breaks_m, genome_2, genome_1, breaks_f).tag = 5;
        

    	//individual undergoes remaining meiosis
        if (individual.tag == 0) {
	    	meiosis_reps = asInteger(spo_spores_per/4)-2
	    	for (rep in 1:reps){
	            breaks = sim.chromosome.drawBreakpoints(individual);
	            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            
	            //save all children produced to vector
	            children = c(child1, child2, child3, child4, child5, child6, child7, child8, child9, child10)
	            children.tag = 0

	            // Mother's fitness affects gametophyte fitness; see survival()
	            if (spo_maternal_effect > 0){
	                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
	            }
	    	}
	    }
    }

    if (individual.tag == 45) {
        //make clones and move to p0 as diploids
        for (i in 1:spo_clones_per)
            p0.addClones(individual).tag = 44;
        
        //individuals also selfs    
        //perform meiosis twice for the selfed diploid

        // 4 gametophytes that producearchegonia+antheridia per meiosis
        breaks_m = sim.chromosome.drawBreakpoints(individual);
        child5 = p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks_m)
        child6 = p0.addRecombinant(NULL, NULL, NULL, genome_1, genome_2, breaks_m)
        child7 = p0.addRecombinant(genome_2, genome_1, breaks_m, NULL, NULL, NULL)

        //4 gametophytes that produce archegonia+antheridia per meiosis
        breaks_f = sim.chromosome.drawBreakpoints(individual);
        child8 = p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks_m)
        child9 = p0.addRecombinant(NULL, NULL, NULL, genome_1, genome_2, breaks_m)
        child10 = p0.addRecombinant(genome_2, genome_1, breaks_m, NULL, NULL, NULL)

        // one gametophyte from each set above produces selfed
        p0.addRecombinant(genome_1, genome_2, breaks_m, genome_2, genome_1, breaks_f).tag = 5;
        

    	//individual undergoes remaining meiosis
    	meiosis_reps = asInteger(spo_spores_per/4)-2
    	for (rep in 1:reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4, child5, child6, child7, child8, child9, child10)
            children.tag = 0

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
    	}   
    }
    
    //sporophyte reproduces normally, adds gametophyte gam-selfing to next gen
    if (individual.tag == 6){ 
    	breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=6;          
            child5 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child6 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child7 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);

        meiosis_reps = asInteger(spo_spores_per/4)-1
    	for (rep in 1:reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4, child5, child6, child7)
            children.tag = 0

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
    	}  
    }
"""

# PARAMETERS
#gam_archegonia_per
# -------------------------
# TAGS
# 0, 1, 2, 4, 5
REPRO_PTER_HOMOSPORE_P0 = """
	//all gametophytes are hermaphroditic
	if (individudal.tag == 0) {
		eggs = gam_archegonia_per
		//for each egg per gametophyte, perform fertilization
		for (rep in 1:am_archegonia_per){
			sperm = p0.sampleIndividuals(1); //all individuals make sperm
            
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                
                // Mother's fitness affects sporophyte fitness; see survival()
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                
                // take out of the mating pool
                sperm.tag = 20;
            }
        }
	}
    
    if (individual.tag == 4) { //add new gametophyte clones to p1 as haploids
        for (i in 1:gam_clones_per)
            p1.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }
    
    if (individual.tag == 44){ //sporophyte clone from last gen moves directly to p1
        //add sporophytic clones
        for (i in 1:spo_clones_per)
            p1.addCloned(individual).tag = 0;
    }
    //move sporophytic selfed into p1
    if (individual.tag == 5)
        p1.addCloned(individual).tag = 0;
    
    if (individual.tag == 6){ //performs gametophytic selfing
        child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL,  NULL);
        child.tag=0;
        // Mother's fitness affects sporophyte fitness; see survival()
        if (gam_maternal_effect > 0)
            child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
    }
"""

LATE_PTER_HOMOSPORE = """
	//tag gametophytes that will self
    p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate));
        clones.tag = 4; //tag clones
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        //tag sporophytes that will clone
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo
        
        //tag sporophytes that will self
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds
        
        //tag sporophytes that will have gametophytes that self
        num_gam_self = rbinom(1, length(p1_size), gam_self_rate);
        gam_selfed_inds = p1.sampleIndividuals(num_gam_self);
        gam_selfed_cloned = gam_selfed_inds[gam_selfed_inds.tag == 44];
        gam_selfed_cloned.tag = 46; //tag selfing and cloning inds
        
        gam_selfed = gam_selfed_inds[gam_selfed_inds.tag == 0];
        gam_selfed.tag = 6; //tag gametophytic selfing 
    }
"""

REPRO_PTER_HOMOSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    if (individual.tag == 0) { //normal sporophyte makes female and male spores
        //chance of creating a merisitc (egg-bearing) gametophyte
        if (runif(1) < gam_female_to_male_ratio){
            meiosis_reps = asInteger(spo_megaspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
                child1.tag = 1;
                
                child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
                child2.tag = 1;
                
                // Mother's fitness affects gametophyte fitness; see survival()
                if (spo_maternal_effect > 0){
                    child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    child2.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
        //else made a male-only gametophyte
        else {
            meiosis_reps = asInteger(spo_microspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
                child1.tag = 2;
                
                child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
                child2.tag = 2;
                
                // Mother's fitness affects gametophyte fitness; see survival()
                if (spo_maternal_effect > 0){
                    child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    child2.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
    }
    
    if (individual.tag == 44) { //save cloned spo to p0
        //make the clonoes
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 44;
        
        //individual can make gametes too
        if (runif(1) < gam_female_to_male_ratio){
            meiosis_reps = asInteger(spo_megaspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
                child1.tag = 1;
                
                child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
                child2.tag = 1;
                
                // Mother's fitness affects gametophyte fitness; see survival()
                if (spo_maternal_effect > 0){
                    child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    child2.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
        
        //else made a male-only gametophyte
        else {
            meiosis_reps = asInteger(spo_microspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
                child1.tag = 2;
                
                child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
                child2.tag = 2;
                
                // Mother's fitness affects gametophyte fitness; see survival()
                if (spo_maternal_effect > 0){
                    child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    child2.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
    }
    
    if (individual.tag == 5) { //sporophytic selfing
        meiosis_reps = asInteger(spo_megaspores_per/2);
        for (rep in 1:meiosis_reps){
            //generate 4 spores (2 rounds of meiosis) with their own breakpoints
            breaks = sim.chromosome.drawBreakpoints(individual); //male
            breaks2 = sim.chromosome.drawBreakpoints(individual); //hermaphrodite
            p0.addRecombinant(NULL, NULL, NULL, g_2, g_1, breaks).tag = 1; // male outcross
            p0.addRecombinant(g_1, g_2, breaks2, NULL, NULL, NULL).tag = 2; //female outcross
            
            p0.addRecombinant(g_1, g_2, breaks, g_2, g_1, breaks2).tag = 5; //add the diploid selfed
        }
        //make the rest of the microspores
        male_meiosis_reps = asInteger(spo_microspores_per/2) - meiosis_reps;
        for (rep in 1:male_meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            for (sperm in 1:gam_sperm_per_microspore){
                p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
                p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 2;
            }
        }
    }

"""

REPRO_PTER_HETEROSPORE_P0 = """
    // females find male gametes to reproduce
    if (individual.tag == 1) { //meristic (egg-bearing) gametophyte
        reproduction_opportunity_count = gam_eggs_per_megaspore;
        for (repro in seqLen(reproduction_opportunity_count)) {
            sperm = p0.sampleIndividuals(1); //all individuals make sperm
            
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                
                // Mother's fitness affects sporophyte fitness; see survival()
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                
                // take out of the mating pool
                sperm.tag = 20;
            }
        }
    }
    
    
    if (individual.tag == 44){ //sporophyte clone from last gen moves directly to p1   
        //add sporophytic clones
        for (i in 1:spo_clones_per)
            p1.addCloned(individual).tag = 0;
    }
    
    //move sporophytic selfed into p1
    if (individual.tag == 5)
        p1.addCloned(individual).tag = 0;
"""

LATE_PTER_HETEROSPORE = """
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        //sporophyte clones
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo
        
        //sporophytic selfing
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds
    }
"""