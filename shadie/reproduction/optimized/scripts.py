#!/usr/bin/env python

"""
String scripts for reproduction blocks
"""

#----------------------------------------
#for reading from population file
READIN_RESCHEDULE = """
finalgen = {sim_time} + sim.generation - 1;
// scheduling the end of the simulation
sim.rescheduleScriptBlock(s0, generations=finalgen);
"""

#----------------------------------------
#early ()
EARLY = """
// diploids (p1) just generated haploid gametophytes
    if (sim.generation % 2 == 0) {{

        // fitness affects gametophyte survival
        megaspores = p0.individuals[p0.individuals.tag==1];
        megaspores.fitnessScaling = (gam_pop_size*(gam_female_to_male_ratio)
                /length(megaspores));
        
        //extra pressure applied to sperm to reduce sim size
        microspores = p0.individuals[p0.individuals.tag==2];
        microspores.fitnessScaling = (microspore_pool/length(microspores));

        //set mutation rate for haploids
        sim.chromosome.setMutationRate(gam_mutation_rate);

        // p0 and p1 survival callbacks
        s1.active = 1;
        s2.active = 0;
        s3.active = 1;
        s4.active = 0;

        // haploids get modified fitness, without dominance
        {activate}
    }}


    // odd generations = gametophytes (p0) just generated sporophytes
    else {{

        // fitness affects sporophytes
        p1.fitnessScaling = spo_pop_size / p1.individualCount;

        //set mutation rate for diploids
        sim.chromosome.setMutationRate(spo_mutation_rate);

        // turn off p0 survival callbacks
        // turn on p1 survival callbacks
        s1.active = 0;
        s2.active = 1;
        s3.active = 0;
        s4.active = 1;

        // diploids get SLiM's standard fitness calculation, with dominance
        {deactivate}
    }}
"""

#-----------------------------------------------
#FITNESS CALLBACKS

#activate/deactivate
ACTIVATE = "{idx}.active = 1;"

DEACTIVATE = "{idx}.active = 0;"

#-----------------------------------------------
#SURVIVAL CALLBACKS
#standard defaults

SPO_DEATH_CHANCE = """
//this code implements random chance of death in sporophytes
    if (runif(1) < spo_random_death_chance)
        return F;
    else
        return NULL;
"""

GAM_DEATH_CHANCE = """
//this code implements random chance of death in gametophytes
    if (runif(1) < gam_random_death_chance)
        return F;
    else
        return NULL;
"""

MATERNAL_EFFECT_P0 = """
    // maternal effect as weighted average
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * gam_maternal_effect) + fitness * (1 - gam_maternal_effect);
        return (draw < corrected_fitness);
    }
"""

MATERNAL_EFFECT_P1 = """
    // maternal effect as weighted average
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * spo_maternal_effect) + fitness * (1 - spo_maternal_effect);
        return (draw < corrected_fitness);
    }
"""

DEATH = """
return F;
"""

#------
#SHADIE SURVIVAL BLOCKS

SURV = """
s1 survival(p1) {{
    return F;
}}

s2 survival(p1) {{
    // this code implements random chance of death in sporophytes
    if (runif(1) < spo_random_death_chance)
        return F;
    
    {p0maternal_effect}
    else
        return NULL;
}}

// even
s3 survival(p0) {{
    //this code implements random chance of death in gametophytes
    if (runif(1) < gam_random_death_chance)
        return F;   
    {p0survival} 
    {p1maternal_effect}
    else
        return NULL;
}}

// odd
s4 survival(p0) {{
    return F;
}}
"""

ANGIO_SURV_P0 = """
    if {//All ovules survive; this is a way of implementing maternal effect
    //if mother died, they would not be produced
        if (individual.tag == 1)
                return T;
    }
"""
#-----------------------------------------------

SUBSTITUTION = """
    // gametophytes have just undergone fitness selection
    if (sim.generation % 2 == 0) {{
        {muts} {late}
"""

#late() **for every mut!
SUB_MUTS = """
        mut{idx} = sim.mutationsOfType({mut});
        freq{idx} = sim.mutationFrequencies(NULL, mut{idx});
        if (any(freq{idx} == 0.5))
            sim.subpopulations.genomes.removeMutations(mut{idx}[freq{idx} == 0.5], T);
"""

#--------------------------
REPRO_BRYO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;
    if (individual.tag == 0) { //normal sporophyte makes female and male spores
        // is it a female gametophyte?
        if (runif(1)<gam_female_to_male_ratio){
            meiosis_reps = asInteger(spo_megaspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                for (egg in 1:gam_eggs_per_megaspore){
                    p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
                    p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
                }
            }
        }
        else { //it is male
            meiosis_reps = asInteger(spo_microspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                for (sperm in 1:gam_sperm_per_microspore){
                    p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 2;
                    p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 2;
                }
            }
        }
    }
    
    if (individual.tag == 4) { //move clones directly to p0
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }
    
    if (individual.tag ==5){
        meiosis_reps = asInteger(spo_megaspores_per/2);
        for (rep in 1:meiosis_reps){
            //generate 4 spores (2 rounds of meiosis) with their own breakpoints
            breaks = sim.chromosome.drawBreakpoints(individual); //male
            breaks2 = sim.chromosome.drawBreakpoints(individual); //female
            p0.addRecombinant(NULL, NULL, NULL, g_2, g_1, breaks).tag = 1; // male outcross
            p0.addRecombinant(g_1, g_2, breaks2, NULL, NULL, NULL).tag = 2; //female outcross
            
            p0.addRecombinant(g_1, g_2, breaks, g_2, g_1, breaks2).tag = 5; //add the diploid selfed
        }
        //make the rest of the sperm
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

REPRO_BRYO_DIO_P0 = """
    // females find male gametes to reproduce
    if (individual.tag == 1) {
        reproduction_opportunity_count = gam_eggs_per_megaspore;
        for (repro in seqLen(reproduction_opportunity_count)) {
            sperm = p0.sampleIndividuals(1, tag=2);
            
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
    
    if (individual.tag == 4) { //add gametophyte clones to p1
        for (i in 1:gam_clones_per){
            p1.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
        }
    }
    
    //move sporophytic selfed into p1
    if (individual.tag == 5) {p1.addCloned(individual).tag=0;}
"""

LATE_BRYO_DIO = """
    p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate));
        clones.tag = 4; //tag clones
    }
	//odd = starts with gam in p0, generates spo into p1
	else {
		p1_size = length(p1.individuals);
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed = p1.sampleIndividuals(number_selfed);
        selfed.tag = 5; //tag sporophytic selfing inds
        
        num_gam_self = rbinom(1, length(p1_size), gam_self_rate);
        gam_selfed = p1.sampleIndividuals(num_gam_self);
        gam_selfed.tag = 6; //tag gametophytic selfing
	}
"""

GAMETOPHYTIC_SELFING = """
    //P1
    if (individual.tag == 6){ //gametophytic selfing
        meiosis_count = asInteger(spo_spores_per/2);
        for (i in 1:meiosis_count)
            breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 6; //one egg will undergo selfing
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag =
            ifelse (runif(1)<gam_female_to_male_ratio, 1, 2); //the other will outcross
    }
    //P0
    if (individual.tag == 6) { //performed gametophytic selfing
    p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL,  NULL).tag=0;
    }
"""

REPRO_BRYO_MONO_P1 = """
    // creation of gametes from sporophytes
    g_1 = genome1;
    g_2 = genome2;

    meiosis_reps = floor(spo_spores_per/2);
    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""

REPRO_BRYO_MONO_P0 = """
    reproduction_opportunity_count = gam_sporophytes_per;

    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= gam_clone_rate)
        {
        reproduction_opportunity_count = reproduction_opportunity_count 
        + (gam_clones_per*gam_sporophytes_per);}

    for (repro in seqLen(reproduction_opportunity_count))
    {
        if (runif(1) <= gam_self_rate)
        {
            // this is selfing using two identical gametes – intragametophytic selfing
            p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
            // intergametophytic selfing might happen below, by chance
        }
        else
        {
            sperm = p0.sampleIndividuals(1);

            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);

            if (gam_maternal_effect > 0) //Mother's fitness affects sporophyte fitness; see survival()
                child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));

        }
    }
"""

LATE_BRYO_MONO = """
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        p1.individuals.tag = 0; //set null tag
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate)); //there is no spo cloning, so remove later
        clones.tag = 4; //tag clones
        
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds
        
        selfed_cloned = selfed_inds[selfed_inds.tag == 4];
        selfed_cloned.tag = 45; //tag selfing and cloning inds
        
        num_gam_self = rbinom(1, length(p1_size), gam_self_rate);
        gam_selfed_inds = p1.sampleIndividuals(num_gam_self);
        
        gam_selfed = gam_selfed_inds[gam_selfed_inds.tag == 0];
        gam_selfed.tag = 6; //tag gametophytic selfing indsp1_size = length(p1.individuals);
        
        gam_selfed_cloned = gam_selfed_inds[gam_selfed_inds.tag == 4];
        gam_selfed_cloned.tag = 46; //tag selfing and cloning inds
    
    }
"""

EARLY1_ANGIO = """
    sim.addSubpop('p1', spo_pop_size); // diploid sporophyte pop
    sim.addSubpop('p0', 0); // haploid gametophyte pop

    fems = spo_female_to_male_ratio*spo_pop_size;
    spo_sex_starts = c(rep(1, asInteger(fems)), 
        rep(0, asInteger(spo_pop_size-fems)));
    p1.individuals.tag = spo_sex_starts;
"""

REPRO_ANGIO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    // individual is female
    if (individual.tag == 1) {

        // determine how many ovules were fertilized, out of the total
        fertilized_ovules = rbinom(1, flower_ovules_per, ovule_fertilization_rate);
        meiosis_reps = floor(fertilized_ovules/2);
        if (runif(1) <= spo_clone_rate)
            meiosis_reps = spo_clones_per*meiosis_reps*2;

        for (rep in 1:meiosis_reps) {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
        }
    }

    // individual is male
    else {
        meiosis_reps = floor(pollen_count/2);
        if (runif(1) <= spo_clone_rate)
            meiosis_reps = spo_clones_per*meiosis_reps*2;
        for (rep in 1:meiosis_reps) {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
        }
    }
"""

LATE_ANGIO_DIO = """
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        p1.individuals.tag = 0; //set null tag
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate)); //there is no spo cloning, so remove later
        clones.tag = 4; //tag clones
        
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds
        
        selfed_cloned = selfed_inds[selfed_inds.tag == 4];
        selfed_cloned.tag = 45; //tag selfing and cloning inds
        
        num_gam_self = rbinom(1, length(p1_size), gam_self_rate);
        gam_selfed_inds = p1.sampleIndividuals(num_gam_self);
        
        gam_selfed = gam_selfed_inds[gam_selfed_inds.tag == 0];
        gam_selfed.tag = 6; //tag gametophytic selfing indsp1_size = length(p1.individuals);
        
        gam_selfed_cloned = gam_selfed_inds[gam_selfed_inds.tag == 4];
        gam_selfed_cloned.tag = 46; //tag selfing and cloning inds
    
    }
"""


REPRO_ANGIO_DIO_P0  = """
    // females find male gametes to reproduce
    if (individual.tag == 1) {
        if (pollen_comp == T) {

            // sperm land on stigma
            pollen_pool = p0.sampleIndividuals(pollen_per_ovule, tag=0);
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

                for (i in range(1:length(pollen_pool)))
                    {attempts = attempts + 1;
                    if (runif(1)<pollen_success_rate)
                        break;
                    }
                idx = attempts-1;    
                target_fitness = sorted_fitness_vector[idx];
                winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
                sperm = winners[0];
            }
            // find a male
            else sperm = p0.sampleIndividuals(1, tag=0);
        }

        else
            // find a male
            sperm = p0.sampleIndividuals(1, tag=0);

        if (sperm.size() == 1) {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
            sperm.tag = 2;

            if (runif(1) <= spo_female_to_male_ratio)
                child.tag = 1;
            else
                child.tag = 0;
        }
    }
"""

REPRO_ANGIO_MONO_P1="""
    g_1 = genome1;
    g_2 = genome2;

    // determine how many ovules were fertilized, out of the total
    fertilized_ovules = rbinom(1, flower_ovules_per, ovule_fertilization_rate);
    meiosis_reps = floor(fertilized_ovules/2);
    if (runif(1) <= spo_clone_rate)
        meiosis_reps = spo_clones_per*meiosis_reps*2;

    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
    }

    meiosis_reps = floor(pollen_count/2);
    if (runif(1) <= spo_clone_rate)
        meiosis_reps = spo_clones_per*meiosis_reps*2;
    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
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
    
    if (individual.tag == 4) //move gam clones directly to p0
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 1;
    
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
    
    if (individual.tag == 45) { //move clones directly to p0
        //make clones
        for (i in 1:spo_clones_per)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 1;
        
        //individuals also selfs    
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
    
    if (individual.tag == 6){ //sporophyte reproduces normally, adds female gametes for gam selfing to next gen
        if (runif(1) < gam_female_to_male_ratio){
            meiosis_reps = asInteger(spo_megaspores_per/2);
            for (rep in 1:meiosis_reps){
                breaks = sim.chromosome.drawBreakpoints(individual);
                child1 = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
                child1.tag = 6;
                
                child2 = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
                child2.tag = 6;
                
                // Mother's fitness affects gametophyte fitness; see survival()
                if (spo_maternal_effect > 0){
                    child1.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    child2.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }
            }
        }
    }
"""

REPRO_PTER_HOMOSPORE_P0 = """
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
    
    
    if (individual.tag == 4) { //add new gametophyte clones to p1
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
    
    if (individual.tag == 45){
        //clones
        for (i in 1:gam_clones_per)
            p1.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
        
        //selfed
        if (individual.tag == 5)
            p1.addCloned(individual).tag = 0;
    }
    
    if (individual.tag == 6){ //performs gametophytic selfing
        child = p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL,  NULL);
        child.tag=0;
        // Mother's fitness affects sporophyte fitness; see survival()
        if (gam_maternal_effect > 0)
            child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
    }
"""

LATE_PTER_HOMOSPORE = """
    p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate));
        clones.tag = 4; //tag clones
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        
        clones = p1.sampleIndividuals(asInteger(p1_size*spo_clone_rate));
        clones.tag = 44; //tag clones - 4 is gam, 44 is spo
        
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed_inds = p1.sampleIndividuals(number_selfed);
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45; //tag selfing and cloning spo inds
        
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5; //tag sporophytic selfing inds
        
        num_gam_self = rbinom(1, length(p1_size), gam_self_rate);
        gam_selfed_inds = p1.sampleIndividuals(num_gam_self);
        gam_selfed_cloned = gam_selfed_inds[gam_selfed_inds.tag == 44];
        gam_selfed_cloned.tag = 46; //tag selfing and cloning inds
        
        gam_selfed = gam_selfed_inds[gam_selfed_inds.tag == 0];
        gam_selfed.tag = 6; //tag gametophytic selfing 
    }
"""

REPRO_PTER_HETEROSPORE_P1 = """
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