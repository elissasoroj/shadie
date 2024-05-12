#!/usr/bin/env python

"""Special case Wright-Fisher specific SLIM script snippets used for
string substitution.
"""

#############################################################
# typical Wright-Fisher with soft selection
# ==PARAMETERS==
# K
# -------------------------
REPRO_WF = """
    // parents are chosen proportional to fitness
    fitness = p2.cachedFitness(NULL);
    parent = sample(p2.individuals, K, replace=T, weights=fitness);
    for (i in seqLen(K))
        p2.addRecombinant(parent.genome1[i], NULL, NULL, NULL, NULL, NULL);

    self.active = 0;
"""

# ==PARAMETERS==
# K
# -------------------------
REPRO_HAPLOID_WF = """
    // parents are chosen randomly. Two haploid genomes
    // come together and immediately produce a haploid child with recombination

    fitness = p2.cachedFitness(NULL);
    parent1 = sample(p2.individuals, K, replace=T);
    parent2 = sample(p2.individuals, K, replace=T);
    for (i in seqLen(K)){
        breaks = sim.chromosome.drawBreakpoints(parent1[i]);
        p2.addRecombinant(parent1.genome1[i], parent2.genome1[i], breaks, NULL, NULL, NULL);
    }

    self.active = 0;
"""

REPRO_HAPLOID_SOFT_WF = """
    // parents are chosen proportional to fitness. Two haploid genomes
    // come together and immediately produce a haploid child with recombination

    fitness = p2.cachedFitness(NULL);
    parent1 = sample(p2.individuals, K, replace=T, weights=fitness);
    parent2 = sample(p2.individuals, K, replace=T, weights=fitness);
    for (i in seqLen(K)){
        breaks = sim.chromosome.drawBreakpoints(parent1[i]);
        p2.addRecombinant(parent1.genome1[i], parent2.genome1[i], breaks, NULL, NULL, NULL);
    }

    self.active = 0;
"""

REPRO_CLONAL_WF = """
    // parents are chosen randomly, produce one offpsring each time.

    fitness = p2.cachedFitness(NULL);
    parent1 = sample(p2.individuals, K, replace=T);

    for (i in seqLen(K)){
        p2.addRecombinant(parent1.genome1[i], NULL, NULL, NULL, NULL, NULL);
    }

    self.active = 0;
"""

REPRO_CLONAL_SOFT_WF = """
    // parents are chosen proportional to fitness, produce one offpsring each time.

    fitness = p2.cachedFitness(NULL);
    parent1 = sample(p2.individuals, K, replace=T, weights=fitness);

    for (i in seqLen(K)){
        p2.addRecombinant(parent1.genome1[i], NULL, NULL, NULL, NULL, NULL);
    }

    self.active = 0;
"""

# PARAMETERS
# GAM_POP_SIZE
# -------------------------
REPRO_ALTGEN_P2 = """
    // parents are chosen randomly
    fitness = p2.cachedFitness(NULL);
    parents = sample(p2.individuals, GAM_POP_SIZE, replace=T);
    for (i in seqLen(GAM_POP_SIZE)){
        breaks = sim.chromosome.drawBreakpoints(parents[i]);
        p1.addRecombinant(parents.genome1[i], parents.genome2[i], breaks, NULL, NULL, NULL);
    }
    self.active = 0;
"""

REPRO_ALTGEN_SOFT_P2 = """
    // parents are chosen proportional to fitness
    fitness = p2.cachedFitness(NULL);
    parents = sample(p2.individuals, GAM_POP_SIZE, replace=T, weights=fitness);
    for (i in seqLen(GAM_POP_SIZE)){
        breaks = sim.chromosome.drawBreakpoints(parents[i]);
        p1.addRecombinant(parents.genome1[i], parents.genome2[i], breaks, NULL, NULL, NULL);
    }
    self.active = 0;
"""

# PARAMETERS
# SPO_POP_SIZE
# -------------------------
REPRO_ALTGEN_P1 = """
    // parents are chosen randomly
    fitness = p1.cachedFitness(NULL);
    parents1 = sample(p1.individuals, SPO_POP_SIZE, replace=T);
    parents2 = sample(p1.individuals, SPO_POP_SIZE, replace=T);
    for (i in seqLen(SPO_POP_SIZE))
        p2.addRecombinant(parents1.genome1[i], NULL, NULL, parents2.genome1[i], NULL, NULL);
    self.active = 0;
"""

REPRO_ALTGEN_SOFT_P1 = """
    // parents are chosen proportional to fitness
    fitness = p1.cachedFitness(NULL);
    parents1 = sample(p1.individuals, SPO_POP_SIZE, replace=T, weights=fitness);
    parents2 = sample(p1.individuals, SPO_POP_SIZE, replace=T, weights=fitness);
    for (i in seqLen(SPO_POP_SIZE))
        p2.addRecombinant(parents1.genome1[i], NULL, NULL, parents2.genome1[i], NULL, NULL);
    self.active = 0;
"""

# ------------------
WF_ALTGEN_EARLY = """
    // diploids (p2) just generated haploid into p1
    if (community.tick % 2 == 0) {
        // fitness affects haploid survival
        p1.fitnessScaling = GAM_POP_SIZE / p1.individualCount;
    }
    // haploids (p1) just generated diploids into p2
    else {
        //fitness affects diploid survival
        p2.fitnessScaling = SPO_POP_SIZE / p2.individualCount;
    }
"""

# ------------
SURV_MORAN = """
    if (individual.age>1)
        return F;
    else
        return NULL;
"""

# ------------
OLD_SURV_WF = """
    if (individual.age>1)
        return F;
    else
        return NULL;
"""
