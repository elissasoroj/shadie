#!/usr/bin/env python

"""
Spermatophyte (just angiosperm currently) string substitutions.
"""

#PARAMETERS
# spo_pop_size
# spo_female_to_male_ratio
# -------------------------
# TAGS
# 0
EARLY1_ANGIO = """
    sim.addSubpop('p1', spo_pop_size);
    sim.addSubpop('p0', 0);

    // tag individuals as male or female.
    fems = spo_female_to_male_ratio * spo_pop_size;
    spo_sex_starts = c(rep(1, asInteger(fems)), 
        rep(2, asInteger(spo_pop_size-fems)));
    p1.individuals.tag = spo_sex_starts;
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
	//sporophytes make hermaphroditic flowers
	if (individual.tag == 1) {
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
	}
	if (individual.tag ==2){
		//calculate chanec a give microsporophyte will be successful
		//total number of ovules in gametophyte population:
		ovule_pop = spo_flowers_per*flower_ovules_per*length(p1.individuals);
		males = length(p1.individuals[p1.individuals.tag==2]);
		successful_pollen = rbinom(1,pollen_comp_stigma_pollen_per*ovule_pop, 1/males);
		pollen_chance = pollen_comp_stigma_pollen_per*ovule_pop/males;
		for (count in 1:successful_pollen){
			
			if (runif(1) < pollen_chance^3){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^2){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^1){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else{
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
			
			}
		}
	}
	
	if (individual.tag == 41) { //save cloned spo to p0
		//make the clones
		for (i in 1:spo_clones_per)
			p0.addCloned(individual).tag = 41;
		
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
	}
	if (individual.tag == 42) { //save cloned spo to p0
		//make the clones
		for (i in 1:spo_clones_per)
			p0.addCloned(individual).tag = 42;
		
		//clone can reproduce normally as well	
		//make microspores (pollen)
		ovule_pop = spo_flowers_per*flower_ovules_per*length(p1.individuals);
		males = length(p1.individuals[p1.individuals.tag==2]);
		successful_pollen = rbinom(1,pollen_comp_stigma_pollen_per*ovule_pop, 1/males);
		pollen_chance = pollen_comp_stigma_pollen_per*ovule_pop/males;
		for (count in 1:successful_pollen){
			
			if (runif(1) < pollen_chance^3){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^2){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^1){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else{
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
			
			}
		}
	}
"""

# PARAMETERS
# pollen_comp_stigma_pollen_per
#spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2
REPRO_ANGIO_DIO_P0  = """
   // females find male gametes to reproduce
	if (individual.tag == 1) {
		if (pollen_comp == T) {
			
			// sperm land on stigma
			pollen_pool = p0.sampleIndividuals(pollen_comp_stigma_pollen_per, tag=2);
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
			p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL).tag= ifelse(runif(1)<spo_female_to_male_ratio, 1, 2);
			sperm.tag = 20; //sperm goes into used pool
		}
	}
	
	//move clones back into p1,reset tag
	if (individual.tag == 41)
		p1.addCloned(individual).tag=1;
	
	if (individual.tag == 42)
		p1.addCloned(individual).tag=2;
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
		maleclones = clones[clones.tag == 2];
		femclones.tag == 41;
		maleclones.tag == 42;
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
		
		//calculate chanec a give microsporophyte will be successful
		//total number of ovules in gametophyte population:
		ovule_pop = spo_flowers_per*flower_ovules_per*length(p1.individuals);
		males = length(p1.individuals[p1.individuals.tag==0]);
		successful_pollen = rbinom(1,pollen_comp_stigma_pollen_per*ovule_pop, 1/males);
		pollen_chance = pollen_comp_stigma_pollen_per*ovule_pop/males;
		for (count in 1:successful_pollen){
			
			if (runif(1) < pollen_chance^3){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^2){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^1){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else{
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
			
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
		ovule_pop = spo_flowers_per*flower_ovules_per*length(p1.individuals);
		males = length(p1.individuals[p1.individuals.tag==0]);
		successful_pollen = rbinom(1,pollen_comp_stigma_pollen_per*ovule_pop, 1/males);
		pollen_chance = pollen_comp_stigma_pollen_per*ovule_pop/males;
		for (count in 1:successful_pollen){
			
			if (runif(1) < pollen_chance^3){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^2){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
			}
			else if (runif(1) < (pollen_chance)^1){
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
			}
			else{
				breaks = sim.chromosome.drawBreakpoints(individual);
				p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
			
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
			
			//only one megaspore will be produced, used for the new selfed sporophyte
			breaks_f = sim.chromosome.drawBreakpoints(individual);
			
			// add the diploid selfed 
			p0.addRecombinant(genome1, genome2, breaks_m, genome2, genome1, breaks_f).tag = 5;
			
			//make extra pollen, if lucky
			//make microspores (pollen)
			ovule_pop = spo_flowers_per*flower_ovules_per*length(p1.individuals);
			males = length(p1.individuals[p1.individuals.tag==0]);
			successful_pollen = rbinom(1,pollen_comp_stigma_pollen_per*ovule_pop, 1/males);
			pollen_chance = pollen_comp_stigma_pollen_per*ovule_pop/males;
			for (count in 1:successful_pollen){
				
				if (runif(1) < pollen_chance^3){
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				}
				else if (runif(1) < (pollen_chance)^2){
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				}
				else if (runif(1) < (pollen_chance)^1){
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				}
				else{
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
				
				}
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
			
			//only one megaspore will be produced, used for the new selfed sporophyte
			breaks_f = sim.chromosome.drawBreakpoints(individual);
			
			// add the diploid selfed 
			p0.addRecombinant(genome1, genome2, breaks_m, genome2, genome1, breaks_f).tag = 5;
			
			//make extra pollen, if lucky
			ovule_pop = spo_flowers_per*flower_ovules_per*length(p1.individuals);
			males = length(p1.individuals[p1.individuals.tag==0]);
			successful_pollen = rbinom(1,pollen_comp_stigma_pollen_per*ovule_pop, 1/males);
			pollen_chance = pollen_comp_stigma_pollen_per*ovule_pop/males;
			for (count in 1:successful_pollen){
				
				if (runif(1) < pollen_chance^3){
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				}
				else if (runif(1) < (pollen_chance)^2){
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
				}
				else if (runif(1) < (pollen_chance)^1){
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag=2;
					p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag=2;
				}
				else{
					breaks = sim.chromosome.drawBreakpoints(individual);
					p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
				
				}
			}
		}
	}
"""

# PARAMETERS
# pollen_comp_stigma_pollen_per
#spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2, 44, 5
REPRO_ANGIO_MONO_P0="""
    // females find male gametes to reproduce
	if (individual.tag == 1) {
		if (pollen_comp == T) {
			
			// sperm land on stigma
			pollen_pool = p0.sampleIndividuals(pollen_comp_stigma_pollen_per, tag=2);
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
		p1.addCloned(individual).tag=0;
	
	//move clones back into p1, reset tag
	if (individual.tag == 5)
		p1.addCloned(individual).tag=0;
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
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds;
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""

