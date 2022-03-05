#!/usr/bin/env python

"""
Special case Wright-Fisher specific SLIM script snippets used for string substitution.
"""

# ==PARAMETERS==
#K
# -------------------------

REPRO_WF = """
    // parents are chosen proportional to fitness
    fitness = p1.cachedFitness(NULL);
    parent = sample(p1.individuals, K, replace=T, weights=fitness);
    for (i in seqLen(K))
        p1.addRecombinant(parent.genome1[i], NULL, NULL, NULL, NULL, NULL);
    
    self.active = 0;
"""

# ==PARAMETERS==
#K
# -------------------------

REPRO_CLONAL_WF = """
    // parents are chosen proportional to fitness. Two haploid genomes 
    come together and immediately produce a haploid child with recombination

    fitness = p1.cachedFitness(NULL);
    parent1 = sample(p1.individuals, K, replace=T, weights=fitness);
    parent2 = sample(p1.individuals, K, replace=T, weights=fitness);
    for (i in seqLen(K)){
        breaks = sim.chromosome.drawBreakpoints(parent1[i]);
        p1.addRecombinant(parent1.genome1[i], parent2.genome1[i], breaks, NULL, NULL, NULL);
    }
    
    self.active = 0;
"""

# PARAMETERS
#GAM_POP_SIZE
# -------------------------
REPRO_ALT_GEN_P1 = """
    // parents are chosen proportional to fitness
    fitness = p1.cachedFitness(NULL);
    parents = sample(p1.individuals, GAM_POP_SIZE, replace=T, weights=fitness);
    for (i in seqLen(GAM_POP_SIZE)){
        breaks = sim.chromosome.drawBreakpoints(parents[i]);
        p0.addRecombinant(parents.genome1[i], parents.genome2[i], breaks, NULL, NULL, NULL);
    }
    self.active = 0;
"""

# PARAMETERS
#GAM_POP_SIZE
# -------------------------
REPRO_ALT_GEN_P0 = """
    // parents are chosen proportional to fitness
    fitness = p0.cachedFitness(NULL);
    parents1 = sample(p0.individuals, SPO_POP_SIZE, replace=T, weights=fitness);
    parents2 = sample(p0.individuals, SPO_POP_SIZE, replace=T, weights=fitness);
    for (i in seqLen(SPO_POP_SIZE))
        p1.addRecombinant(parents1.genome1[i], NULL, NULL, parents2.genome1[i], NULL, NULL);
    self.active = 0;
"""

#------------
OLD_SURV_WF = """
    if (individual.age>1) 
        return F;
    else
        return NULL;
"""