#!/usr/bin/env python

"""
Bryophyte specific SLIM script snippets used for string substitution.
"""

# PARAMETERS
# gam_female_to_male_ratio
# spo_spores_per
#gam_archegonia_per
#gam_antheridia_per
# -------------------------
# TAGS
# 0, 1, 2, 4, 5
REPRO_BRYO_DIO_P1 = """
    // normal sporophyte: make both female and male spores
    if (individual.tag == 0) { //focus on one sporophyte

        // each meiotic division produces 4 spores, two female and two male
        meiosis_reps = asInteger(spo_spores_per/4);
        for (rep in 1:meiosis_reps) {

            // sample a meiosis crossover position.
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create 2 megaspores and 2 microspores (all recombinant for now)
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 2;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 2;
        }
    }

    // cloned gametophytes from last gen: copy and move clones directly to p0,
    //these cloned gametophyte becomes a 'normal' gametophyte in p0, tag=4 removed
    if (individual.tag == 41) 
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 1;

    if (individual.tag == 42) 
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 2;

    // selfing sporophyte: make male, female, and self-recombinant spores
    if (individual.tag == 5) {
        // perform rounds of two meiotic divisions
        breaks1 = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 1;
        microspore1 = p0.addRecombinant(genome1, genome2, breaks1, NULL, NULL, NULL);
        microspore2 = p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL);

        breaks2 = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(genome1, genome2, breaks2, NULL, NULL, NULL).tag = 1;
        p0.addRecombinant(genome2, genome1, breaks2, NULL, NULL, NULL).tag = 1;
        microspore_3 = p0.addRecombinant(genome2, genome1, breaks2, NULL, NULL, NULL);

        // add the diploid selfed 
        p0.addRecombinant(genome1, genome2, breaks1, genome2, genome1, breaks2).tag = 5;

        //see note below**
        microspores = c(microspore1, microspore2, microspore3);
        microspores.tag = 2;

        // perform any additional meiosis rounds for the rest of the spores
        // each meiotic division produces 4 spores, two female and two male
        meiosis_reps = asInteger(spo_spores_per/4)-8;
        for (rep in 1:meiosis_reps) {

            // sample a meiosis crossover position.
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create 2 megaspores and 2 microspores (all recombinant for now)
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 1;
            microspore1 = p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL);
            microspore2 = p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL);

            //**each microspore will produce a male gametophyte that may
            //develop many antheridia, each of which produces thousands of sperm.
            //For simplicity and computation time we do not model these sperm individuallally,
            //but this part of the code may be modified to implement parameters such as 
            //chance of male gametophyte developing antheridium, etc. 
            microspores = c(microspore1, microspore2);
            microspores.tag = 2;
        }
    }
"""

# Parameters:
# gam_archegonia_per (we could use "megaspore_archegonia_per" or "gam_eggs_per")
# gam_clones_per
# -------------------------
# TAGS:
# 0, 1, 4, 5, 20
REPRO_BRYO_DIO_P0 = """
    // P0 at this point contains microspores and megaspores
    // normal gametophyte: if female, fertilize with recombinant male
    if (individual.tag == 1) {

    	// each megaspore makes N archegonia, each of which contains one egg
        eggs = gam_archegonia_per;
        for (egg in 1:eggs) {

            // sample a sperm from p0; each microspore creates one sperm
            sperm = p0.sampleIndividuals(1, tag=2);
            //note, the sperm will NOT be removed from the pool, because it is actually
            //a pool of microspores, which can produce thousands of clonal sperm
            
            //if we run out of sperm, skip zygote formation
            if (sperm.size() == 1) {

                // create zygote and set tag to 0
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag = 0;

                // store maternal fitness to apply maternal effects later (in survival).
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
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
# spo_self_rate
# ------------------
# TAGS
# 4, 5
LATE_BRYO_DIO = """
        p0_size = length(p0.individuals);
        clones = p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate));
        femclones = clones[clones.tag ==1];
        femclones.tag = 41; //tag female clones;
        maleclones = clones[clones.tag ==2];
        maleclones.tag = 42; //tag male clones;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed = p1.sampleIndividuals(number_selfed);
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""


## PARAMETERS
# spo_spores_per
# ------------------
# TAGS
# 0, 4, 5
REPRO_BRYO_MONO_P1 = """
    // normal sporophyte: make both female and male spores
    if (individual.tag == 0) { //focus on one sporophyte

        // each meiotic division produces 4 spores
        meiosis_reps = asInteger(spo_spores_per/4);
        for (rep in 1:meiosis_reps) {

            // sample a meiosis crossover position.
            breaks = sim.chromosome.drawBreakpoints(individual);

            // create 2 megaspores and 2 microspores (all recombinant for now)
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(genome1, genome2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(genome2, genome1, breaks, NULL, NULL, NULL).tag = 0;
        }
    }

    // cloned gametophytes from last gen: copy and move clones directly to p0,
    //these cloned gametophyte becomes a 'normal' gametophyte in p0, tag=4 removed
    if (individual.tag == 4) 
        p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 0;

    // selfing sporophyte: make selfed each meiotic division
    if (individual.tag == 5) {
        meiosis_reps = asInteger(spo_spores_per/4);
        for (rep in 1:meiosis_reps)
        breaks1 = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome1, genome2, breaks1, NULL, NULL, NULL).tag = 0;p0.addRecombinant(genome2, genome1, breaks1, NULL, NULL, NULL).tag = 0;

        breaks2 = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(genome1, genome2, breaks2, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome1, genome2, breaks2, NULL, NULL, NULL).tag = 0;
        p0.addRecombinant(genome2, genome1, breaks2, NULL, NULL, NULL).tag = 0;

        // add the diploid selfed 
        p0.addRecombinant(genome1, genome2, breaks1, genome2, genome1, breaks2).tag = 5;

        }
"""

## PARAMETERS
# 
# ------------------
# TAGS
# 0, 4, 5
REPRO_BRYO_MONO_P0 = """
    // p0 at this point contains microspores and megaspores
    // normal gametophyte: makes both archegonia and antheridia
    if (individual.tag == 0) {

        // each megaspore makes N archegonia, each of which contains one egg
        eggs = gam_archegonia_per;
        for (egg in 1:eggs) {

            // sample a sperm from p0; each microspore creates many sperm
            sperm = p0.sampleIndividuals(1);
            //note, the sperm will NOT be removed from the pool, because it is actually
            //a pool of microspores, which can produce thousands of clonal sperm
            
            if (sperm.size() == 1) {

                // create zygote and set tag to 0
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag = 0;

                // store maternal fitness to apply maternal effects later (in survival).
                if (gam_maternal_effect > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
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
    if (individual.tag == 5) 
        p1.addCloned(individual).tag = 0;
"""

LATE_BRYO_MONO = """
//odd = starts with gam in p0, generates spo into p1
    p0_size = length(p0.individuals);
    //tag gametophytic clones
    p0.sampleIndividuals(asInteger(p0_size*gam_clone_rate)).tag=4;
    }
    //odd = starts with gam in p0, generates spo into p1
    else {
        p1_size = length(p1.individuals);
        number_selfed = rbinom(1, length(p1_size), spo_self_rate);
        selfed = p1.sampleIndividuals(number_selfed);
        selfed.tag = 5; //tag sporophytic selfing inds;
    }
"""
