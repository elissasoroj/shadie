#!/usr/bin/env python

"""Spermatophyte (just angiosperm currently) string substitutions.
"""

# shadie DEFINITIONS
DEFS_ANGIO_MONO = """
// model: monecious angiosperm
// p1 = haploid population
// p2 = diploid population
// 1 = gametophyte (1N) tag
// 2 = gametophyte clones (1N) tag ----------- none in this model
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag

// L0: Female = T, Male = F
"""

DEFS_ANGIO_DIO = """
// model: dioecious angiosperm
// p0 = pollen
// p1 = haploid female population
// p2 = diploid population
// 1 = gametophyte (1N) tag
// 2 = gametophyte clones (1N) tag ----------- none in this model
// 3 = sporophyte (2N) tag
// 4 = sporophyte clones (2N) tag

// L0: Female = T, Male = F
"""

# -------------

P1_ANGIO_FITNESS_SCALE = """
    // random death also occurs to implement GAM_CEILING
    if (p0.individualCount > GAM_CEILING) {
        to_kill = p0.individualCount-GAM_CEILING;
        random_death = sample(p0.individuals, to_kill);
        sim.killIndividuals(p0.individuals[random_death.index]);
    }

    // fitness affects female gametophyte survival
    p1.fitnessScaling = GAM_POP_SIZE / p1.individualCount;
"""

P2_ANGIO_FITNESS_SCALE = """
    //kill off pollen
    sim.killIndividuals(p0.individuals);

    // fitness is scaled relative to number of inds in p2
    p2.fitnessScaling = SPO_POP_SIZE / p2.individualCount;
"""

EGG_FITNESS_AFFECTS_VIABILITY = """
    // determine how many ovules were fertilized, out of the total
    ind_fitness = p2.cachedFitness(ind.index);
    max_fitness = max(p2.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;
    
    meiosis_reps = rbinom(1, SPO_ARCHEGONIA_PER, ind_fitness_scaled);
"""

NO_EGG_FITNESS = "meiosis_reps = SPO_ARCHEGONIA_PER;"

FITNESS_AFFECTS_POLLEN_NUM = """
    // determine how many pollen grains produced
    ind_fitness = p2.cachedFitness(ind.index);
    max_fitness = max(p2.cachedFitness(NULL));
    ind_fitness_scaled = ind_fitness/max_fitness;
    
    pollen = rbinom(1, SPO_POLLEN_PER, ind_fitness_scaled);
"""

CONSTANT_POLLEN_NUM = "pollen = SPO_POLLEN_PER;"


POLLEN_COMPETITION = """
    // sperm land on stigma
    pollen_pool = sample(males, STIGMA_POLLEN_PER);
    for (pollen in pollen_pool) {
        // store fitness value
        pollen.setValue("fitness", p0.cachedFitness(pollen.index));
    }

    if (length(pollen_pool)>0) {
        //sort pollens by fitness
        fitness_vector = pollen_pool.getValue("fitness");
        sorted_fitness_vector = sort(fitness_vector, ascending=F);
    
        //calculate how many pollens attempt to fertilize
        attempts = 0;
        for (i in range(1:length(pollen_pool))) {
            attempts = attempts + 1;
            if (runif(1)<POLLEN_SUCCESS_RATE)
                break;
        }
        idx = attempts-1;
        target_fitness = sorted_fitness_vector[idx];
        winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
        sperm = winners[0];
    }
"""

NO_POLLEN_COMPETITION = """
        //no pollen competition
        sperm = sample(males, 1);
"""

FIRST1_ANGIO_MONO = """
    sim.addSubpop('p2', SPO_POP_SIZE);
    sim.addSubpop('p1', 0);
    sim.addSubpop('p0', 0); // pollen
    p2.individuals.tag = 3;
"""

FIRST1_ANGIO_DIO = """
    sim.addSubpop('p2', SPO_POP_SIZE);
    sim.addSubpop('p1', 0); //eggs
    sim.addSubpop('p0', 0); // pollen
    p2.individuals.tag = 3;
    p2.individuals.setValue('maternal_fitness', 1.0);
    p2.individuals.tagL0 = (runif(p2.individualCount) < SPO_FEMALE_TO_MALE_RATIO);
"""


EARLY_ANGIO_DIO = """
    // random death also occurs to implement GAM_CEILING
    if (p0.individualCount > GAM_CEILING) {
        to_kill = GAM_CEILING - GAM_POP_SIZE;
        death_chance = to_kill/p0.individualCount;
        random_death = sample(c(F,T), p0.individualCount, T, c(1-death_chance, death_chance));
        sim.killIndividuals(p0.individuals[random_death]);
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
# 1, 2, 41, 42
REPRO_ANGIO_DIO_P2 = """
    ind = individual;
    
    // clonal individual get added to the p1 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < SPO_CLONE_RATE) {{
        for (i in seqLen(SPO_CLONES_PER)) {{
            child = p2.addRecombinant(ind.genome1, NULL, NULL, ind.genome2, NULL, NULL, parent1=ind);
            child.tag = 4; // sporophyte clone
            child.tagL0 = ind.tagL0;
            
            //sporophyte maternal effect not applicable for clones = neutral
         if (SPO_MATERNAL_EFFECT > 0)
            child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }}
    }}
    
    //only females will make eggs
    if (individual.tagL0) {{
    {egg_production}
    
        //one egg per archegonia (fertilized ovule)
        for (rep in seqLen(meiosis_reps)){{
            breaks = sim.chromosome.drawBreakpoints(individual);
            egg = p1.addRecombinant(ind.genome1, ind.genome2, breaks, NULL, NULL, NULL, parent1=ind);
            egg.tag = 1;
            egg.tagL0 = T;
            if (SPO_MATERNAL_EFFECT > 0)
                egg.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }}
    }}

    //males will make pollen
    else {{
        {pollen_production}
        meiosis_reps = asInteger(floor(pollen/4));
        for (rep in seqLen(meiosis_reps)){{
            breaks1 = sim.chromosome.drawBreakpoints(ind);
            breaks2 = sim.chromosome.drawBreakpoints(ind);
        
            // create four meiotic products
            child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL, parent1=ind);
            child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL, parent1=ind);
            child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL, parent1=ind);
            child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL, parent1=ind);
            children = c(child1, child2, child3, child4);
            children.tag = 1;       
            children.tagL0 = F;
            if (SPO_MATERNAL_EFFECT > 0)
                children.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }}
    }}
"""

# PARAMETERS
# stigma_pollen_per
# spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2
REPRO_ANGIO_DIO_P1 = """
    ind = individual;

    //reproduction scripts run only females
    if (individual.tagL0) {{
    
    // get all males that could fertilize an egg of this female
    males = p0.individuals;
    
    // Each egg is outcrossed in this model
    {pollen_competition}
    child = p2.addRecombinant(ind.genome1, NULL, NULL, sperm.genome1, NULL, NULL, parent1=ind, parent2=sperm);
    child.tag = 3;
    child.tagL0 = ifelse(runif(1) > SPO_FEMALE_TO_MALE_RATIO, T, F);
    }}
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
REPRO_ANGIO_MONO_P2 = """
    ind = individual;
    
    // clonal individual get added to the p1 pool for next round.
    // NOTE: this doesn't allow clones to reproduce this round.
    if (runif(1) < SPO_CLONE_RATE) {{
        for (i in seqLen(SPO_CLONES_PER)) {{
            child = p2.addRecombinant(ind.genome1, NULL, NULL, ind.genome2, NULL, NULL, parent1=ind);
            child.tag = 4; // sporophyte clone
            
            //sporophyte maternal effect not applicable for clones = neutral
         if (SPO_MATERNAL_EFFECT > 0)
            child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
        }}
    }}
    
    {egg_production}
    
    //one egg per archegonia (fertilized ovule)
    for (rep in seqLen(meiosis_reps)){{
        breaks = sim.chromosome.drawBreakpoints(individual);
        egg = p1.addRecombinant(ind.genome1, ind.genome2, breaks, NULL, NULL, NULL, parent1=ind);
        egg.tag = 1;
        egg.tagL0 = T;
    }}
    
    {pollen_production}

    meiosis_reps = asInteger(floor(pollen/4));
    for (rep in seqLen(meiosis_reps))
    {{
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        
        // create four meiotic products
        child1 = p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL, parent1=ind);
        child2 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL, parent1=ind);
        child3 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL, parent1=ind);
        child4 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL, parent1=ind);
        children = c(child1, child2, child3, child4);
        children.tag = 1;       
        children.tagL0 = F;
    }}
    if (SPO_MATERNAL_EFFECT > 0)
        children.setValue("maternal_fitness", subpop.cachedFitness(ind.index));
"""

# PARAMETERS
# stigma_pollen_per
# spo_female_to_male_ratio
# -------------------------
# TAGS
# 1, 2, 44, 5

REPRO_ANGIO_MONO_P1 = """
    ind = individual;
    
    //reproduction scripts run only females
    if (ind.tagL0) {{
    
    // get all males that could fertilize an egg of this female
    males = p0.individuals[p0.individuals.tagL0==F];
    
    // if selfing is possible then get all sibling males
   if (SPO_SELF_RATE_PER_EGG > 0)
    //shared parent count for sibs is > 0 (1 or 2)
    siblings = males[individual.sharedParentCount(males)!=0];
    
    //find non-related individuals
    outcross_males = males[individual.sharedParentCount(males)==0];
    
    // weighted sampling: each egg is spo-selfed, or outcrossed.
    //gametophytic selfing cannot occur in this model
    mode = sample(
    x=c(1, 2, 3),
      size=1,
      weights=c(0, SPO_SELF_RATE_PER_EGG, 1 - (SPO_SELF_RATE_PER_EGG))
      );    
      
      if (mode == 2) {{
        if (siblings.size() > 0) {{
            sibling = sample(siblings, 1);
            child = p2.addRecombinant(ind.genome1, NULL, NULL,
                sibling.genome1, NULL, NULL, parent1 = individual, parent2 = sibling);
            child.tag = 3; //sporophyte tag
                }}
            }}
        if (mode == 3) {{
            {pollen_competition}
            child = p2.addRecombinant(ind.genome1, NULL, NULL, sperm.genome1, NULL, NULL, parent1=ind, parent2=sperm);
            child.tag = 3;
        }}
    }}
"""
