#!/usr/bin/env python

"""Pteridophyte-specific scripts.
"""

EARLY_PTER_HETERO = """
    // reproduction(p1) just finished creating sperm and eggs in p0.
    // this generation will form sporophytes and apply selection on them.
    if (sim.generation % 2 == 0) {{

        // fitness is scaled relative to non-males in p0
        num_males = length(p0.individuals[p0.individuals.tag == 2]);
        p0.fitnessScaling = (GAM_POP_SIZE / (p0.individualCount - num_males));
        // print(format('p0.fitnessScaling=%.2f', p0.fitnessScaling));

        // set mutation rate to gametophyte rate
        sim.chromosome.setMutationRate(GAM_MUTATION_RATE);

        // turn on p0 and off p1 survival callbacks
        s1.active = 1;
        s2.active = 0;
        s3.active = 1;
        s4.active = 0;

        // haploids reproduce, diploids don't
        s5.active = 0;
        s6.active = 1;

        // activate/deactivate mutations that apply to haploids
        // and turns off any dominance effects
        {mutations_activate}
    }}


    // reproduction(p0) just finished creating diploid zygotes
    // this generation will produce gametophytes and apply selection on them.
    else {{

        // fitness is scaled relative to number of inds in p1
        p1.fitnessScaling = SPO_POP_SIZE / p1.individualCount;
        // print(format('p1.fitnessScaling=%.2f', p1.fitnessScaling));

        // set mutation rate to sporophyte rate
        sim.chromosome.setMutationRate(SPO_MUTATION_RATE);

        // turn off p0 and on p1 survival callbacks
        s1.active = 0;
        s2.active = 1;
        s3.active = 0;
        s4.active = 1;

        // diploids reproduce, haploids don't
        s5.active = 1;
        s6.active = 0;

        // activate/deactivate mutations that apply fitness to diploids
        // and turns on any dominance effects
        {mutations_deactivate}
    }}
"""


# This is called at the end of every generation to remove fixed mutations,
# tag ...
LATE_PTER_HETERO = """\
    // gametophytes have just undergone fitness selection
    if (sim.generation % 2 == 0) {{

        // get the total number of genomes (p1 is empty right now)
        p0_size = p0.individualCount;

        // for each MutationType check all muts for fixation
        {checking_each_mut_for_fixation}

        // tag gametophytes that will clone
        clones = p0.sampleIndividuals(asInteger(p0_size * GAM_CLONE_RATE));
        clones.tag = 4;
    }}

    // sporophytes have just undergone fitness selection
    else {{
        p1_size = length(p1.individuals);
        // print(format('p1_size=%i', p1_size));

        // tag sporophytes that will clone (44)
        clones = p1.sampleIndividuals(asInteger(p1_size * SPO_CLONE_RATE));
        clones.tag = 44;

        // tag sporophytes that will clone and self (45)
        selfed_inds = p1.sampleIndividuals(asInteger(p1_size * SPO_SELF_RATE_PER_IND));
        selfed_cloned = selfed_inds[selfed_inds.tag == 44];
        selfed_cloned.tag = 45;

        // tag sporophytes that will self (5)
        selfed = selfed_inds[selfed_inds.tag == 0];
        selfed.tag = 5;
    }}
"""


# FUNCTIONS SPECIFIC TO PTER_HETERO:
# make_microspores(ind, reps)
# make_eggs(ind, reps)
# sporophyte_selfs(ind)
FUNCTIONS_PTER_HETERO = """
// p0 = haploid population
// p1 = diploid population

// 0 = hermaphrodite
// 1 = female
// 2 = male
// 4 = gametophyte clone
// 44 = sporophyte clone
// 45 = sporophyte cloned and selfed
// 5 = sporophytic selfed
// 6 = gametophyte selfed


// shadie func(ind, reps) to generate 4 p0 tag=2 (male/microspores) in heterosporous ferns.
// Each ind here represents a 2N microsporangium, and each replicate is a different
// mitotic 2N sporocyte that it produces. Each sporocyte undergoes meiosis to
// form 4 1N recombinant gametes (microspores) that are added to p0.
function (void)make_microspores(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){

        // 4 microspores per meiosis rep
        breaks1 = sim.chromosome.drawBreakpoints(individual);
        breaks2 = sim.chromosome.drawBreakpoints(individual);
        child1 = p0.addRecombinant(individual.genome1, individual.genome2, breaks1, NULL, NULL, NULL);
        child2 = p0.addRecombinant(individual.genome2, individual.genome1, breaks1, NULL, NULL, NULL);
        child3 = p0.addRecombinant(individual.genome1, individual.genome2, breaks2, NULL, NULL, NULL);
        child4 = p0.addRecombinant(individual.genome2, individual.genome1, breaks2, NULL, NULL, NULL);

        // create vector of all recombined children genomes
        children = c(child1, child2, child3, child4);
        children.tag = 2;

        // set fitness so that mother affects gametophyte fitness; see survival()
        if (SPO_MATERNAL_EFFECT > 0)
            children.setValue("maternal_fitness", individual.subpopulation.cachedFitness(individual.index));
    }
}


// shadie func(ind, reps) to generate 1 p0 tag=1 (female/egg) in heterosporous fern
function (void)make_eggs(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){
        breaks = sim.chromosome.drawBreakpoints(individual);
        child1 = p0.addRecombinant(individual.genome1, individual.genome2, breaks, NULL, NULL, NULL);
        child1.tag = 1;

        // set fitness so that mother affects gametophyte fitness; see survival()
        if (SPO_MATERNAL_EFFECT > 0)
            child1.setValue("maternal_fitness", individual.subpopulation.cachedFitness(individual.index));
    }
}


// shadie func(ind) to sporophytic self in heterosporous fern
// generates p1 tag=5 (self) selfed sporophytes
// generates p0 tag=0 (herm) gametophytes, and p0 tag=1 or tag=2 (eggs, sperm)
function (void)sporophyte_selfs(object<Individual>$ ind) {

    // count the number of eggs that are selfed, only the unselfed ones
    // will be added to the p0 pool.
    selfed_eggs = 0;

    // iterate over megaspores per megasporangia
    for (i in 1:SPO_MEGASPORES_PER) {

        // heterosporous ferns have 1-few archegonia (and thus egg) per megasporangium
        for (i in 1:GAM_ARCHEGONIA_PER) {

            // randomly sample whether this egg was self fertilized
            if (runif(1) < SPO_SELF_RATE_PER_EGG) {
                selfed_eggs = selfed_eggs + 1;

                // create 4 self sperm, use one of them to create new p1 zygote,
                // put other 3 sperm into the gamete pool.
                breaks1 = sim.chromosome.drawBreakpoints(individual);
                breaks2 = sim.chromosome.drawBreakpoints(individual);
                breaks_f = sim.chromosome.drawBreakpoints(individual);
                p1.addRecombinant(ind.genome1, ind.genome2, breaks1, ind.genome1, ind.genome2, breaks_f).tag = 5;

                // TODO: if these are extra sperm why are they not tag=2?
                // 4 gametophytes that produce antheridia per meiosis
                child1 = p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL);
                child2 = p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL);
                child3 = p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL);

                // save all children produced to vector
                children = c(child1, child2, child3);
                children.tag = 0;

                // mother's fitness affects gametophyte fitness; see survival()
                if (SPO_MATERNAL_EFFECT > 0)
                    children.setValue("maternal_fitness", ind.subpopulation.cachedFitness(individual.index));
            }
            else
                make_eggs(individual, 1);
        }
    }

    // gametophyte is male: make microspores minus the ones already used in selfing.
    make_microspores(individual, SPO_MICROSPORES_PER - selfed_eggs);
}
"""


REPRO_PTER_HETEROSPORE_P1 = """
    // sporophyte is hermaphroditic non-clonal and non-selfing
    if (individual.tag == 0) {
        make_eggs(individual, SPO_MEGASPORES_PER);
        make_microspores(individual, max(4, SPO_MICROSPORES_PER));
    }

    // sporophyte creates clones for next generation
    // and also creates eggs & sperms this generation
    if (individual.tag == 44) {
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        make_eggs(individual, SPO_MEGASPORES_PER);
        make_microspores(individual, max(4, SPO_MICROSPORES_PER));
    }

    // sporophytic selfed only
    if (individual.tag == 5)
        sporophyte_selfs(individual);

    // sporophytic cloned and selfed
    if (individual.tag == 45) {
        for (i in 1:SPO_CLONES_PER)
            p1.addCloned(individual).tag = 44;
        sporophyte_selfs(individual);
    }
"""


# p0 gametophytes are either female (1) or male (2)
REPRO_PTER_HETEROSPORE_P0 = """
    // spore is female (megaspore)
    if (individual.tag == 1) {

        // iterate over each egg to fertilize.
        for (rep in 1:GAM_ARCHEGONIA_PER) {

            // sample a microspore
            // NOTE: each microspore gives rise to a male gametophyte, which will
            // produce many antheridia, giving rise to thousands of clonal sperm
            // because of this, sperm is not removed from the mating pool when used
            sperm = p0.sampleIndividuals(1, tag=2);

            // add new p1 sporophyte tag=0 (hermaph)
            if (sperm.size() == 1) {
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                child.tag=0;
            }
        }
    }

    // if gametophyte clone (tag=4), add p0 gametophyte clone
    else if (individual.tag == 4) {
        for (i in 1:GAM_CLONES_PER)
            p0.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL).tag = 4;
    }
"""


# SURVIVAL functions are similar for all clades except in the specific
# tags that are referenced for whether the taxon is a spo or gam.
SURVIVAL_PTER_HETERO = """
// remove p1 individuals during even generations
s1 survival(p1) {
    if ((individual.tag == 44) | (individual.tag == 5) | (individual.tag == 45)) {
        individual.tag = 0;
        return T;
    }
    else
        return F;
}

// remove p1s by random chance of death, or by maternal effects to fitness
s2 survival(p1) {
    if (runif(1) < SPO_RANDOM_DEATH_CHANCE)
        return F;

    // gametophyte mother fitness affects sporophyte survival
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * GAM_MATERNAL_EFFECT) + fitness * (1 - GAM_MATERNAL_EFFECT);
        return (runif(1) < corrected_fitness);
    }
    return NULL;
}

// remove p0s by random chance of death and apply maternal effects to fitness
s3 survival(p0) {
    //this code implements random chance of death in gametophytes
    if (runif(1) < GAM_RANDOM_DEATH_CHANCE)
        return F;

    // sporophyte mother fitness affects gametophyte survival
    maternal_effect = individual.getValue("maternal_fitness");
    if (!isNULL(maternal_effect)) {
        corrected_fitness = (maternal_effect * SPO_MATERNAL_EFFECT) + fitness * (1 - SPO_MATERNAL_EFFECT);
        return (runif(1) < corrected_fitness);
    }
    return NULL;
}

// remove p0 individuals during odd generations
s4 survival(p0) {
    if ((individual.tag == 4) | (individual.tag == 6)) {
        individual.tag = 0;
        return T;
    }
    return F;
}
"""
