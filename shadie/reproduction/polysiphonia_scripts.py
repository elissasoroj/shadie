#!/usr/bin/env python

"""Polysiphonia specific SLIM script snippets used for string substitution.

...
"""

PSIPH_FITNESS_SCALE = """p1.fitnessScaling = (GAM_POP_SIZE*(1/GAM_FEMALE_TO_MALE_RATIO))/
(GAM_ARCHEGONIA_PER*p1.individualCount);"""

DEFS_PSIPH = """
//Polysiphonia - triphasic
// p0 = mixed population (caroposporophytes)
// p1 = haploid population (gametophytes)
// p2 = diploid population (tetrasporophytes)

// 1 = gametophyte (1N) tag
// 2 = gametophyte clones (1N) tag
// 3 = tetrasporophyte (2N) tag
// 4 = tetrasporophyte clones (2N) tag;
// 5 = carposporophyte (2N) tag

// L0: Female = True, Male = False
"""

# ==PARAMETERS==
#TSPO_CLONE_RATE
#TSPO_CLONES_PER
#GAM_FEMALE_TO_MALE_RATIO
# -------------------------
# TAGS

REPRO_PSIPH_P2 = """
    ind = individual;
    
    // clonal individual get added to the p0 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < TSPO_CLONE_RATE) {{
        for (i in 1:TSPO_CLONES_PER) {{
            child = p2.addRecombinant(individual.genome1, NULL, NULL, individual.genome2, NULL, NULL);
            child.tag = 4;
            child.setValue("parentid", individual.tag);
        }}
    }}

    {tetraspore_determination}

    // each spore produces its own recombinant breakpoints
    for (rep in seqLen(tetraspores)) {{
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);

        // create four meiotic products

        child1 = p1.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL, parent1 = ind);
        child2 = p1.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL, parent1 = ind);
        child3 = p1.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL, parent1 = ind);
        child4 = p1.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL, parent1 = ind);
        children = c(child1, child2, child3, child4);
        children.tag = 1; //gametophyte tag
        children.tagL0 = (runif(4) < GAM_FEMALE_TO_MALE_RATIO); 
    }}

"""

# PARAMETERS
#GAM_CLONE_RATE
#GAM_CLONES_PER
#GAM_MATERNAL_EFFECT
#SPO_SELF_RATE_PER_EGG
#SPO_SELF_RATE
# -------------------------
# TAGS
# 2, >1M-
REPRO_PSIPH_P1 = """
    // clonal individual get added to the p1 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < GAM_CLONE_RATE) {{
        for (i in seqLen(GAM_CLONES_PER)) {{
            child = p1.addRecombinant(individual.genome1, NULL, NULL,
            NULL, NULL, NULL, parent1 = individual); //only 1 parent recorded
            child.tag = 2; //marked as clone
            child.tagL0 = individual.tagL0; //inherits parent sex
            //gametophyte maternal effect not applicable for clones = neutral
            if (GAM_MATERNAL_EFFECT > 0)
                child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }}
    }}

    // iterate over each egg to find a mate (self, sib, or outcross)
    // NOTE: each gametophyte gives rise to antheridia that produce thousands of
    // clonal sperm. Because of this, sperm is not removed from the mating pool when used.

    //Reproduction scripts run only on female gametophytes
    //Females: L0 = T
    if (individual.tagL0) {{

        // get all males that could fertilize an egg of this female
        males = p1.individuals[p1.individuals.tagL0 == F];
        
        // if selfing is possible then get all sibling males
        if (TSPO_SELF_RATE_PER_EGG > 0)
            siblings = males[individual.sharedParentCount(males)!=0];

        {carpogonia_determination}

        // iterate over each reproductive opportunity (carpogonium) on this female.
        for (rep in seqLen(carpogonia)) {{

            // weighted sampling: each egg gets gam-selfed, spo-selfed, or outcrossed.
            // Note: shadie enforces that GAM_SELF_RATE + TSPO_SELF_RATE !> 1
            mode = sample(
                x=c(1, 2, 3),
                size=1,
                weights=c(0.0, TSPO_SELF_RATE_PER_EGG, (1 - TSPO_SELF_RATE_PER_EGG))
                );

            // intra-gametophytic selfed
            if (mode == 1) {{
                child = p0.addRecombinant(individual.genome1, NULL, NULL,
                individual.genome1, NULL, NULL, parent1 = individual, parent2 = individual);
                child.tag = 5; //carposporophyte tag
                //gametophyte maternal effect on new sporophyte
                if (GAM_MATERNAL_EFFECT > 0)
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            }}

            // inter-gametophytic selfed individual (same sporo parent)
            // only occurs IF a sibling gametophyte is still alive.
            else if (mode == 2) {{
                if (siblings.size() > 0) {{
                    sibling = sample(siblings, 1);
                    child = p0.addRecombinant(individual.genome1, NULL, NULL,
                    sibling.genome1, NULL, NULL, parent1 = individual, parent2 = sibling);
                    child.tag = 5; //carposporophyte tag
                    //gametophyte maternal effect on new sporophyte
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                }}
            }}
            // outcrossing individual samples any other p1, and checks that it is outcrossing
            // only occurs if a non-sib gametophyte is still alive.

            else {{
                outcross_sperms = males[individual.sharedParentCount(males)==0];
                if (! isNULL(outcross_sperms)) {{
                    {sperm_sampling}
                    child = p0.addRecombinant(individual.genome1, NULL, NULL,
                    sperm.genome1, NULL, NULL, parent1 = individual, parent2=sperm);
                    child.tag = 5; //carposporophyte tag
                    //gametophyte maternal effect on new sporophyte
                    if (GAM_MATERNAL_EFFECT > 0)
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));


                }}
            }}

        }}
    }}
"""
CONSTANT_CARPOGONIA = "carpogonia = GAM_MAX_CARPOGONIA_PER;"

MATERNAL_FITNESS_AFFECTS_CARPOGONIA_NUM = """
//calculate carpogonia based on maternal fitness
    ind_fitness = p1.cachedFitness(ind.index);
    max_fitness = max(p1.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;
        
    carp_vector = sample(c(0,1), GAM_MAX_CARPOGONIA_PER, replace = T, weights = c((1-ind_fitness_scaled), ind_fitness_scaled));
    carpogonia = sum(carp_vector);
"""

REPRO_PSIPH_P0 = """
    ind = individual;

    {carpospore_determination}

    if (CSPO_RECOMBINATION) {{
        //mitotic recombination based on Bulankova et al. (2021)
        mitoses = floor(carpospores/2);
        for (rep in seqLen(mitoses)) {{
            breaks1 = sim.chromosome.drawBreakpoints(ind);

            child1 = p1.addRecombinant(ind.genome1, ind.genome2, breaks1, ind.genome1, NULL, NULL, parent1 = ind);
            child2 = p1.addRecombinant(ind.genome2, ind.genome1, breaks1, ind.genome2, NULL, NULL, parent1 = ind);
            children = c(child1, child2);
            child.tag = 3; //tetrasporophyte tag
        }}
    }}
    
    else {{
        for (rep in seqLen(carpospores)) {{
            child = p2.addRecombinant(individual.genome1, NULL, NULL, 
                individual.genome2, NULL, NULL, parent1 = individual);
            child.tag = 3; //tetrasporophyte tag
        }}
    }}
"""

CONSTANT_CSPORES = "carpospores = CSPO_MAX_SPORES_PER;"

CARPOGONIUM_FITNESS_AFFECTS_CARPOSPORE_NUM = """
spore_vector = sample(c(0,1), CSPO__MAX_SPORES_PER, replace = T, 
    weights = c((1-ind_fitness_scaled), ind_fitness_scaled));
carpospores = sum(spore_vector);
"""

