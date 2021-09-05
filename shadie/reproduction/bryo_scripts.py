#!/usr/bin/env python

"""
Bryophyte specific SLIM script snippets used for string substitution.
"""

# PARAMETERS
# gam_female_to_male_ratio
# spo_megaspores_per
# gam_eggs_per_megaspore
# spo_microspores_per
# -------------------------
# TAGS
# 0, 1, 2, 4, 5
REPRO_BRYO_DIO_P1 = """
    // normal sporophyte: make both female and male spores
    if (individual.tag == 0) {

        // if it is a female gametophyte
        if (runif(1) < gam_female_to_male_ratio) {

            // for each megaspore in the female gametophyte
            meiosis_reps = asInteger(spo_megaspores_per / 2);
            for (rep in 1:meiosis_reps){

                // sample a meiosis crossover position.
                breaks = sim.chromosome.drawBreakpoints(individual);

                // create N recombinant eggs per, set tag to 1.
                for (egg in 1:gam_eggs_per_megaspore){
                    p0.addRecombinant(genome_1, genome_2, breaks, NULL, NULL, NULL).tag = 1;
                    p0.addRecombinant(genome_2, genome_1, breaks, NULL, NULL, NULL).tag = 1;
                }
            }
        }

        // if it is a male gametophyte
        else {

            // for each microspore in the male gametophyte
            meiosis_reps = asInteger(spo_microspores_per / 2);
            for (rep in 1:meiosis_reps) {

                // sample a meiosis crossover position.
                breaks = sim.chromosome.drawBreakpoints(individual);

                // create N recombinant sperm, set tag to 2.
                for (sperm in 1:gam_sperm_per_microspore){
                    p0.addRecombinant(genome_1, genome_2, breaks, NULL, NULL, NULL).tag = 2;
                    p0.addRecombinant(genome_2, genome_1, breaks, NULL, NULL, NULL).tag = 2;
                }
            }
        }

    // clonal sporophyte: move clones directly to p0, sets tag to 4.
    if (individual.tag == 4) {
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }

    // selfed sporophyte: make male, female, and self-recombinant spores
    if (individual.tag == 5) {

        // for each megaspore in the female gametophyte
        meiosis_reps = asInteger(spo_megaspores_per / 2);
        for (rep in 1:meiosis_reps) {

            // sample meiosis crossover position to generate 4 spores (2 rounds of meiosis)
            // male outcross
            breaks_m = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(NULL, NULL, NULL, genome_2, genome_1, breaks).tag = 1;

            // female outcross
            breaks_f = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(genome_1, genome_2, breaks2, NULL, NULL, NULL).tag = 2;

            // add the diploid selfed
            p0.addRecombinant(genome_1, genome_2, breaks_m, genome_2, genome_1, breaks_f).tag = 5;
        }

        // for each microspore in the male gametophyte
        male_meiosis_reps = asInteger(spo_microspores_per / 2) - meiosis_reps;
        for (rep in 1:male_meiosis_reps){

            // sample a meiosis crossover position
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create N recombinant sperm, set tag to 2.
            for (sperm in 1:gam_sperm_per_microspore) {
                p0.addRecombinant(genome_1, genome_2, breaks, NULL, NULL, NULL).tag = 2;
                p0.addRecombinant(genome_2, genome_1, breaks, NULL, NULL, NULL).tag = 2;
            }
        }
    }
"""

# PARAMETERS:
# gam_eggs_per_megaspore (was gam_sporophytes_per)
# gam_clones_per
# -------------------------
# TAGS:
# 0, 1, 4, 5, 20
REPRO_BRYO_DIO_P0 = """
    // P0 at this point contains microsporangia and megasporangia
    // gametophytic normal: if female, sample new recombinant male and female gametes
    if (individual.tag == 1) {

    	// for each egg (reproductive opportunity) in this indiv.
        reproduction_opportunity_count = gam_eggs_per_megaspore;
        for (repro in seqLen(reproduction_opportunity_count)) {

            // sample a sperm from p0
            sperm = p0.sampleIndividuals(1, tag=2);
            if (sperm.size() == 1) {

                // create zygote and set tag to 0
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag = 0;

                // store maternal fitness to apply maternal effects later (in survival).
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));

                // remove this sperm from the p0 mating pool
                sperm.tag = 20;
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
    if (individual.tag == 5) {
        p1.addCloned(individual).tag = 0;
    }
"""

# PARAMETERS
# gam_clone_rate
#
# ------------------
# TAGS
# 4, 6
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
# gam_maternal_effect (was Maternal_weight)
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
            if (gam_maternal_effect > 0)
                child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }
    }
"""
