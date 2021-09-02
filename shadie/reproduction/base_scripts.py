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
        p0.fitnessScaling = (gam_ne / p0.individualCount);

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
        p1.fitnessScaling = spo_ne / p1.individualCount;

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

MATERNAL_EFFECT = """
    // maternal effect as weighted average
    maternal_effect = individual.getValue("maternal_fitness");

    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * Maternal_weight) + fitness * (1 - Maternal_weight);
        return (draw < corrected_fitness);
    }

    return NULL;
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
    else
        return NULL;
    
    {p1maternal_effect}
}}

// even
s3 survival(p0) {{
    //this code implements random chance of death in gametophytes
    if (runif(1) < gam_random_death_chance)
        return F;
    else
        {p0survival}
    
    {p0maternal_effect}
}}

// odd
s4 survival(p0) {{
    return F;
}}
"""

ANGIO_SURV_P0 = """
    {//All ovules survive; this is a way of implementing maternal effect
    //if mother died, they would not be produced
        if (individual.tag == 1)
                return T;
        else
            return NULL;
    }
"""
#-----------------------------------------------

SUBSTITUTION = """
    // gametophytes have just undergone fitness selection
    if (sim.generation % 2 == 0) {{
        {inner}
    }}
"""

#late() **for every mut!
SUB_INNER = """
        mut{idx} = sim.mutationsOfType({mut});
        freq{idx} = sim.mutationFrequencies(NULL, mut{idx});
        if (any(freq{idx} == 0.5))
            sim.subpopulations.genomes.removeMutations(mut{idx}[freq{idx} == 0.5], T);
"""

#--------------------------
REPRO_BRYO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    meiosis_reps = floor(spores_per_spo/2);
    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 
            ifelse (runif(1)<gam_female_to_male_ratio, 1, 0);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 
            ifelse (runif(1)<gam_female_to_male_ratio, 1, 0);
    }
"""

REPRO_BRYO_DIO_P0 = """
    // females find male gametes to reproduce
    if (individual.tag == 1) {
        reproduction_opportunity_count = spo_per_gam;

        // clones give the focal individual extra opportunities to reproduce
        if (runif(1) <= gam_clone_rate)
        {
            reproduction_opportunity_count = reproduction_opportunity_count 
                + (gam_clone_number*spo_per_gam);
        }

        for (repro in seqLen(reproduction_opportunity_count)) {
            if (runif(1) <= gam_self_rate) {
                // this is selfing using two identical gametes – intragametophytic selfing
                // intergametophytic (sporophytic) selfing might happen below, by chance
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
            }
            else {
                // find a male!
                sperm = p0.sampleIndividuals(1, tag=0);

                if (sperm.size() == 1) {
                    child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                    // Mother's fitness affects sporophyte fitness; see survival()
                    if (gam_maternal_effect > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));

                    // take out of the mating pool
                    sperm.tag = 2;
                }
            }
        }
    }
"""

REPRO_BRYO_MONO_P1 = """
    // creation of gametes from sporophytes
    g_1 = genome1;
    g_2 = genome2;

    meiosis_reps = floor(spores_per_spo/2);
    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""

REPRO_BRYO_MONO_P0 = """
    reproduction_opportunity_count = spo_per_gam;

    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= gam_clone_rate)
        {
        reproduction_opportunity_count = reproduction_opportunity_count 
        + (gam_clone_number*spo_per_gam);}

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

            if (Maternal_weight > 0) //Mother's fitness affects sporophyte fitness; see survival()
                child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));

        }
    }
"""

EARLY1_ANGIO = """
    sim.addSubpop('p1', spo_ne); // diploid sporophyte pop
    sim.addSubpop('p0', 0); // haploid gametophyte pop

    dsex_starts = c(rep(1, asInteger(spo_female_to_male_ratio*spo_ne)), 
        rep(0, asInteger((1-spo_female_to_male_ratio)*spo_ne)));
    p1.individuals.tag = dsex_starts;
"""

REPRO_ANGIO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    // individual is female
    if (individual.tag == 1) {

        // determine how many ovules were fertilized, out of the total
        fertilized_ovules = rbinom(1, ovule_count, ovule_fertilization_rate);
        meiosis_reps = floor(fertilize_ovules/2);
        if (runif(1) <= spo_clone_rate)
            meiosis_reps = spo_clone_number*meiosis_reps*2;

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
            meiosis_reps = spo_clone_number*meiosis_reps*2;
        for (rep in 1:meiosis_reps) {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
        }
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
                sorted_pool = sort(pollen_pool, ascending=F);
                                
                //calculate how many pollens attempt to fertilize
                attempts = 0;

                for (i in range(1:length(pollen_pool)))
                    {attempts = attempts + 1;
                    if (runif(1)<pollen_success_rate)
                        break;
                    }
                winner = attempts-1;    
                sperm = sorted_pool[winner];
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
    fertilized_ovules = rbinom(1, ovule_count, ovule_fertilization_rate);
    meiosis_reps = floor(fertilized_ovules/2);
    if (runif(1) <= spo_clone_rate)
        meiosis_reps = spo_clone_number*meiosis_reps*2;

    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
    }

    meiosis_reps = floor(pollen_count/2);
    if (runif(1) <= spo_clone_rate)
        meiosis_reps = spo_clone_number*meiosis_reps*2;
    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
    }
"""

REPRO_PTER_HOMOSPORE_P0 = """
    // chance of making meristic (egg-bearing) gametophyte
    if (runif(1) < gam_female_to_male_ratio) {
        reproduction_opportunity_count = 1;

        // clones give the focal individual extra opportunities to reproduce
        if (runif(1) <= gam_clone_rate)
            reproduction_opportunity_count = reproduction_opportunity_count + gam_clone_number;

        for (repro in seqLen(reproduction_opportunity_count)) {
            if (runif(1) <= gam_self_rate)
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
                // this is selfing using two identical gametes – intragametophytic selfing
                // intergametophytic (sporophytic) selfing might happen below, by chance

            else {
                sperm = p0.sampleIndividuals(1);

                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);

                //Mother's fitness affects sporophyte fitness; see survival()
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));

            }
        }
    }
"""

REPRO_PTER_HOMOSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    meiosis_reps = floor(spores_per_spo/2);
    reproduction_opportunity_count = 1;

    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= spo_clone_rate)
        meiosis_reps = spo_clone_number*meiosis_reps*2;

    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""

REPRO_PTER_HETEROSPORE_P0 = """
    g_1 = genome1;
    g_2 = genome2;
    if (individual.tag == 1) // reproduction callbacks for megaspores only
    {   reproduction_opportunity_count = 1;
            
        for (repro in seqLen(reproduction_opportunity_count))
        {   sperm = p0.sampleIndividuals(1, tag=0); // find a male!
                
            if (sperm.size() == 1)
            {// intergametophytic/sporophytic selfing might happen by chance
                child = p1.addRecombinant(individual.genome1, NULL, NULL,
                        sperm.genome1, NULL, NULL);
            }
        }
    }
"""

REPRO_PTER_HETEROSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    meiosis_reps = floor(spores_per_spo/2);
    reproduction_opportunity_count = 1;
    
    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= spo_clone_rate)
        meiosis_reps = spo_clone_number*meiosis_reps*2;

    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
	}
	
	if (runif(1) <= gam_female_to_male_ratio)
		child.tag = 1;
    else
        child.tag = 0;
    
    //Mother's fitness affects gametophyte fitness; see survival()
    if (gam_maternal_effect > 0)
        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
"""