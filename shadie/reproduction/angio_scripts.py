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


#
# gam_megagametophyte_number (was ovule_count)
# gam_megagametophyte_fertilization_rate   (was ovule_fertilization_rate)
# spo_clone_rate
# spo_clone_number
# 
REPRO_ANGIO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    // individual is female
    if (individual.tag == 1) {

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


# 
# pollen_comp
# pollen_per_ovule
# pollen_success_rate
# spo_female_to_male_ratio
#
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


# 
# ovule_count
# ovule_fertilization_rate
# spo_clone_rate
# spo_clone_number
# pollen_count
# 
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


#
#
ANGIO_SURV_P0 = """
    {
    // All ovules survive; this is a way of implementing maternal effect
    // if mother died, they would not be produced
        if (individual.tag == 1)
            return T;
        else
            return NULL;
    }
"""
