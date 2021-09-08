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
# 0, 1, 2, 4, 44, 5, 45, 6, 46
REPRO_PTER_HOMOSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    //normal sporophyte makes spores
    if (individual.tag == 0) {
    	meiosis_reps = asInteger(spo_spores_per/4);
    	for (rep in 1:meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4);
            children.tag = 0;

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
	    	meiosis_reps = asInteger(spo_spores_per/4);
	    	for (rep in 1:meiosis_reps){
	            breaks = sim.chromosome.drawBreakpoints(individual);
	            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            
	            children = c(child1, child2, child3, child4);
	            children.tag = 0;

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
        child1 = p0.addRecombinant(NULL, NULL, NULL, genome2, genome1, breaks_m);
        child2 = p0.addRecombinant(NULL, NULL, NULL, genome1, genome2, breaks_m);
        child3 = p0.addRecombinant(genome2, genome1, breaks_m, NULL, NULL, NULL);

        //4 gametophytes that produce archegonia+antheridia per meiosis
        breaks_f = sim.chromosome.drawBreakpoints(individual);
        child4 = p0.addRecombinant(NULL, NULL, NULL, genome2, genome1, breaks_m);
        child5 = p0.addRecombinant(NULL, NULL, NULL, genome1, genome2, breaks_m);
        child6 = p0.addRecombinant(genome2, genome1, breaks_m, NULL, NULL, NULL);

        // one gametophyte from each set above produces selfed
        p0.addRecombinant(genome1, genome2, breaks_m, genome2, genome1, breaks_f).tag = 5;

        //save all children produced to vector
        children = c(child1, child2, child3, child4, child5, child6);
        children.tag = 0;

        // Mother's fitness affects gametophyte fitness; see survival()
        if (spo_maternal_effect > 0){
            children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }
        

    	//individual undergoes remaining meiosis
        if (individual.tag == 0) {
	    	meiosis_reps = asInteger(spo_spores_per/4)-2;
	    	for (rep in 1:meiosis_reps){
	            breaks = sim.chromosome.drawBreakpoints(individual);
	            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
	            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	            
	            //save all children produced to vector
	            children = c(child1, child2, child3, child4);
	            children.tag = 0;

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
            p0.addClone(individual).tag = 44;
        
        //individuals also selfs    
        //perform meiosis twice for the selfed diploid

        // 4 gametophytes that producearchegonia+antheridia per meiosis
        breaks_m = sim.chromosome.drawBreakpoints(individual);
        child1 = p0.addRecombinant(NULL, NULL, NULL, genome2, genome1, breaks_m);
        child2 = p0.addRecombinant(NULL, NULL, NULL, genome1, genome2, breaks_m);
        child3 = p0.addRecombinant(genome2, genome1, breaks_m, NULL, NULL, NULL);

        //4 gametophytes that produce archegonia+antheridia per meiosis
        breaks_f = sim.chromosome.drawBreakpoints(individual);
        child4 = p0.addRecombinant(NULL, NULL, NULL, genome2, genome1, breaks_m);
        child5 = p0.addRecombinant(NULL, NULL, NULL, genome1, genome2, breaks_m);
        child6 = p0.addRecombinant(genome2, genome1, breaks_m, NULL, NULL, NULL);

        // one gametophyte from each set above produces selfed
        p0.addRecombinant(genome1, genome2, breaks_m, genome2, genome1, breaks_f).tag = 5;
        
        children = c(child1, child2, child3, child4, child5, child6);
        children.tag = 0;

        // Mother's fitness affects gametophyte fitness; see survival()
        if (spo_maternal_effect > 0){
            children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }

    	//individual undergoes remaining meiosis
    	meiosis_reps = asInteger(spo_spores_per/4)-2;
    	for (rep in 1:meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4);
            children.tag = 0;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
    	}   
    }
    
    //sporophyte reproduces normally, adds one gametophyte gam-selfing to next gen
    if (individual.tag == 6){ 
    	breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=6;          
        child5 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
        child6 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
        child7 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);

        meiosis_reps = asInteger(spo_spores_per/4)-1;
    	for (rep in 1:meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4, child5, child6, child7);
            children.tag = 0;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
    	}  
    }

    //sporophyte reproduces normally, adds gametophyte gam-selfing to next gen,
    //and also clones
    if (individual.tag == 46){ 
    	//make clones and move to p0 as diploids
        for (i in 1:spo_clones_per)
            p0.addClone(individual).tag = 44;

    	breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=6;          
        child5 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
        child6 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
        child7 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);

        meiosis_reps = asInteger(spo_spores_per/4)-1;
    	for (rep in 1:meiosis_reps){
            breaks = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            child3 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);          
            child4 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
            
            children = c(child1, child2, child3, child4, child5, child6, child7);
            children.tag = 0;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
    	}  
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
		eggs = gam_archegonia_per;
		//for each egg per gametophyte, perform fertilization
		for (rep in 1:gam_archegonia_per){
			sperm = p0.sampleIndividuals(1); //all individuals make sperm;
            //NOTE: each gametophyte makes many antheridia, each of which 
            //produec thousands of clonal sperm. For this reason, sperm is 
            //not removed from the mating pool once used
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                
                // Mother's fitness affects sporophyte fitness; see survival()
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }
	}
    
    if (individual.tag == 4) { //add new gametophyte clones to p1 as haploids
        for (i in 1:gam_clones_per)
            p1.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }
    
    if (individual.tag == 44) //sporophyte clone from last gen moves directly to p1
        p1.addCloned(individual).tag = 0;

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
	//tag gametophytes that will self
    p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate));
        clones.tag = 4; //tag clones;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        //tag sporophytes that will clone
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo;
        
        //tag sporophytes that will self
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds;
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds;
        
        //tag sporophytes that will have gametophytes that self
        num_gam_self = rbinom(1, length(p1_size), gam_self_rate);
        gam_selfed_inds = p1.sampleIndividuals(num_gam_self);
        gam_selfed_cloned = gam_selfed_inds[gam_selfed_inds.tag == 44];
        gam_selfed_cloned.tag = 46; //tag selfing and cloning inds;
        
        gam_selfed = gam_selfed_inds[gam_selfed_inds.tag == 0];
        gam_selfed.tag = 6; //tag gametophytic selfing;
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
    g_1 = genome1;
    g_2 = genome2;

    //normal sporophyte makes female and male spores
    if (individual.tag == 0) { 
        //chance of creating megaspores
        strobilus_female_ratio = rs_megasporangia_per/rs_microsporangia_per;
        if (runif(1) < strobilus_female_ratio){
            meiosis_reps = asInteger(megasporangia_megaspores_per*rs_megasporangia_per);
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
        //else make microspores
        else {
        	//4 microspores per meiosis rep
            meiosis_reps = asInteger((microsporangia_microspores_per*strobilus_microspores_per)/4);
            for (rep in 1:meiosis_reps){
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
    }
    
    if (individual.tag == 44) {
        //make sporophytic clones and save to p0
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 44;
        
        //individual can make gametes too
        //chance of creating megaspores
        strobilus_female_ratio = rs_megasporangia_per/rs_microsporangia_per;
        if (runif(1) < strobilus_female_ratio){
        	megaspores = asInteger(megasporangia_megaspores_per*rs_megasporangia_per)
            meiosis_reps = megaspores;
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
                child1.tag = 1;
                
                // Mother's fitness affects gametophyte fitness; see survival()
                if (spo_maternal_effect > 0){
                    child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    child2.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
        //else make microspores
        else {
        	//4 microspores per meiosis rep
            meiosis_reps = asInteger((microsporangia_microspores_per*strobilus_microspores_per)/4);
            for (rep in 1:meiosis_reps){
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
    }
    
    if (individual.tag == 5) { //sporophytic selfing
        // for each megaspore generated by the sporophyte, perform meiosis twice
        megaspores = asInteger(megasporangia_megaspores_per*rs_megasporangia_per);
        meiosis_reps = megaspores;
        for (rep in 1:meiosis_reps) {

            // sample meiosis crossover position to generate 4 microspores
            // male outcross
            breaks_m = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(NULL, NULL, NULL, genome2, genome1, breaks_m).tag = 2;
            child2 = p0.addRecombinant(NULL, NULL, NULL, genome1, genome2, breaks_m).tag = 2;
            child3 = p0.addRecombinant(genome2, genome1, breaks_m, NULL, NULL, NULL).tag = 2;

            //only one megaspore will be produced, used for the new selfed sporophyte
            breaks_f = sim.chromosome.drawBreakpoints(individual);
            
            // add the diploid selfed 
            p0.addRecombinant(genome1, genome2, breaks_m, genome2, genome1, breaks_f).tag = 5;

            //maternal effects
            children = c(child1, child2, child3);
            children.tag = 2;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }
        }

        // perform any additional meiosis rounds for male
        microspores = microsporangia_microspores_per*rs_microsporangia_per;
        male_meiosis_reps = asInteger(microspores/4) - (meiosis_reps*4);
        for (rep in 1:male_meiosis_reps){

            // sample a meiosis crossover position
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create 4 recombinant sperm, set tag to 2.
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

   	if (individual.tag == 45) { //sporophytic selfing and cloning
        //make sporophytic clones and save to p0
        for (i in 1:spo_clones_per)
            p0.addCloned(individual).tag = 44;

        // for each megaspore generated by the sporophyte, perform meiosis twice
        megaspores = asInteger(megasporangia_megaspores_per*rs_megasporangia_per);
        meiosis_reps = megaspores;
        for (rep in 1:meiosis_reps) {

            // sample meiosis crossover position to generate 4 microspores
            // male outcross
            breaks_m = sim.chromosome.drawBreakpoints(individual);
            child1 = p0.addRecombinant(NULL, NULL, NULL, genome2, genome1, breaks_m).tag = 2;
            child2 = p0.addRecombinant(NULL, NULL, NULL, genome1, genome2, breaks_m).tag = 2;
            child3 = p0.addRecombinant(genome2, genome1, breaks_m, NULL, NULL, NULL).tag = 2;

            //only one megaspore will be produced, used for the new selfed sporophyte
            breaks_f = sim.chromosome.drawBreakpoints(individual);
            
            // add the diploid selfed 
            p0.addRecombinant(genome1, genome2, breaks_m, genome2, genome1, breaks_f).tag = 5;

            //maternal effects
            children = c(child1, child2, child3, child4, child5, child6, child7);
            children.tag = 2;

            // Mother's fitness affects gametophyte fitness; see survival()
            if (spo_maternal_effect > 0){
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }

        }

        // perform any additional meiosis rounds for male
        microspores = microsporangia_microspores_per*rs_microsporangia_per;
        male_meiosis_reps = asInteger(microspores/4) - meiosis_reps;
        for (rep in 1:male_meiosis_reps){

            // sample a meiosis crossover position
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create 4 recombinant sperm, set tag to 2.
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
#gam_archegonia per (female gam = mature megaspore)
# -------------------------
# TAGS
#0, 1, 2, 44, 5, 20
REPRO_PTER_HETEROSPORE_P0 = """
    // find the megaspores
    if (individual.tag == 1) {
        eggs = gam_archegonia_per;
        //fertilize each egg
        for (rep in 1:eggs)) {
            sperm = p0.sampleIndividuals(1, tag==2); //find a microspore
            //NOTE: each microspore gives rise to a male gametophyte, which will
            //produce many antheridia, giving rise to thousands of clonal sperm
            //because of this, sperm is not removed from the mating pool when used
            
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
                
            }
        }
    }
    
    if (individual.tag == 44) //sporophyte clone from last gen moves directly to p1   
        p1.addCloned(individual).tag = 0;

    //move sporophytic selfed into p1
    if (individual.tag == 5)
        p1.addCloned(individual).tag = 0;
"""

# PARAMETERS
#spo_clone_rate
#spo_self_rate
# -------------------------
# TAGS
# 44, 45, 5
LATE_PTER_HETEROSPORE = """
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        //sporophyte clones
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo;
        
        //sporophytic selfing
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds;
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""