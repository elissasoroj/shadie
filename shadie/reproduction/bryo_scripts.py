#!/usr/bin/env python

"""
Bryophyte specific SLIM script snippets used for string substitution.
"""

#
# spo_spores_per 
# gam_female_to_male_ratio
#
REPRO_BRYO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;

    // divided by 2 because ...
    meiosis_reps = floor(spo_spores_per/2);
    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 
            ifelse (runif(1)<gam_female_to_male_ratio, 1, 0);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 
            ifelse (runif(1)<gam_female_to_male_ratio, 1, 0);
    }
"""

#
# gam_sporophytes_per
# gam_clone_rate
# gam_clone_number
# gam_self_rate
# 
REPRO_BRYO_DIO_P0 = """
    // females find male gametes to reproduce
    if (individual.tag == 1) {

    	// add additional reproductive opportunities to indiv if clonal
        reproduction_opportunity_count = gam_sporophytes_per;
        if (runif(1) <= gam_clone_rate) {
            reproduction_opportunity_count = reproduction_opportunity_count 
                + (gam_clone_number * gam_sporophytes_per);
        }

        // iterate over each opportunity to create recombinant gametophytes
        for (repro in seqLen(reproduction_opportunity_count)) {

			// intragametophytic selfing (two identical gametes)
            if (runif(1) <= gam_self_rate) {
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
            }

            // intergametophytic (sporophytic) selfing.
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


# 
# spo_spores_per
# 
REPRO_BRYO_MONO_P1 = """
    // creation of gametes from sporophytes
    g_1 = genome1;
    g_2 = genome2;

    meiosis_reps = floor(spo_spores_per/2);
    for (rep in 1:meiosis_reps) {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""


# 
# gam_sporophytes_per (was spo_per_gam)
# gam_clone_rate
# gam_clone_number
# gam_self_rate
# spo_maternal_effect (was Maternal_weight; should it be gam_maternal_effect?)
# 
REPRO_BRYO_MONO_P0 = """
    reproduction_opportunity_count = gam_sporophytes_per;

    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= gam_clone_rate) {
        reproduction_opportunity_count = reproduction_opportunity_count 
        + (gam_clone_number * gam_sporophytes_per);
    }

    for (repro in seqLen(reproduction_opportunity_count)) {

	    // this is selfing using two identical gametes – intragametophytic selfing
        if (runif(1) <= gam_self_rate) {          
            p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
        }

        // intergametophytic selfing might happen below, by chance
        else {
            sperm = p0.sampleIndividuals(1);
            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
            
            // Mother's fitness affects sporophyte fitness; see survival()
            if (spo_maternal_effect > 0)
                child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }
    }
"""
