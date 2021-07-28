#!/usr/bin/env python

"""
String scripts for reproduction blocks
"""
#----------------------------------------
#early ()
EARLY = """
if (sim.generation % 2 == 0) //diploids (p1) just generated haploid gametophytes
{{
    //fitness affects gametophyte survival
    p0.fitnessScaling = (hK / p0.individualCount);

    //p0 and p1 survival callbacks
    s1.active = 1;
    s2.active = 0;
    s3.active = 1;
    s4.active = 0;

    // haploids get modified fitness, without dominance
    {activate}
}}
else //odd generations = gametophytes (p0) just generated sporophytes
{{
    p1.fitnessScaling = dK / p1.individualCount; //fitness affects sporophytes

    //turn off p0 survival callbacks
    //turn on p1 survival callbacks
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

#callback
FIT = """
return 1 + mut.selectionCoeff; //gametophytes have no dominance effects
"""

#activate/deactivate
ACTIVATE = "{idx}.active = 1;"

DEACTIVATE = "{idx}.active = 0;"

#-----------------------------------------------
#SURVIVAL CALLBACKS
#standard defaults

DEATH_CHANCE = """
//this code implements random death chance
    if (runif(1) < Death_chance)
        return F;
    else
        return NULL;
"""

MATERNAL_EFFECT = """
    // maternal effect as weighted average
    maternal_effect = individual.getValue("maternal_fitness");
    
    if (!isNULL(maternal_effect))
    {
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

SURV ="""
s1 survival(p1)
{{
    return F;
}}

s2 survival(p1)
{{
    //this code implements random death chance
    if (runif(1) < Death_chance)
        return F;
    else
        return NULL;
    
    {maternal_effect}
}}

s3 survival(p0) //even
{{
    //this code implements random death chance
    if (runif(1) < Death_chance)
        return F;
    else
        return NULL;
}}

s4 survival(p0) //odd
{{
    return F;
}}
"""

#-----------------------------------------------
#late() **for every mut!
SUB_INNER = """
mut{idx} = sim.mutationsOfType({mut});
    freq{idx} = sim.mutationFrequencies(NULL, mut{idx});
    if (any(freq{idx} == 0.5))
        sim.subpopulations.genomes.removeMutations(mut{idx}[freq{idx} == 0.5], T)
"""

SUBSTITUTION = """
    if (sim.generation % 2 == 0) //gametophytes have just undergone fitness selection
    {{
        {inner}
    }}
"""
#--------------------------
REPRO_BRYO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;
    
    meiosis_reps = floor(Spore_num/2);
    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = ifelse (runif(1)<FtoM, 1, 0);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = ifelse (runif(1)<FtoM, 1, 0);
    }
"""

REPRO_BRYO_DIO_P0 = """
    if (individual.tag == 1)    // females find male gametes to reproduce
    {
        reproduction_opportunity_count = 1;
        
        // clones give the focal individual extra opportunities to reproduce
        if (runif(1) <= Clone_rate)
            reproduction_opportunity_count = reproduction_opportunity_count + 1;
        
        for (repro in seqLen(reproduction_opportunity_count))
        {
            if (runif(1) <= Self_rate)
            {
                // this is selfing using two identical gametes – intragametophytic selfing
                // intergametophytic selfing might happen below, by chance
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
            }
            else
            {
                sperm = p0.sampleIndividuals(1, tag=0); // find a male!
                
                if (sperm.size() == 1)
                {
                    child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                    
                    if (Maternal_weight > 0) //Mother's fitness affects sporophyte fitness; see survival()
                        child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
                    
                    sperm.tag = 2;  // take out of the mating pool
                }
            }
        }
    }
"""

REPRO_BRYO_MONO_P1 = """
    // creation of gametes from sporophytes
    g_1 = genome1;
    g_2 = genome2;
    
    meiosis_reps = floor(Spore_num/2);
    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""

REPRO_BRYO_MONO_P0 = """
    reproduction_opportunity_count = 1;
    
    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= Clone_rate)
        reproduction_opportunity_count = reproduction_opportunity_count + Clone_num;
    
    for (repro in seqLen(reproduction_opportunity_count))
    {
        if (runif(1) <= Self_rate)
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
    sim.addSubpop('p1', dK); //diploid sporophyte pop
    sim.addSubpop('p0', 0); //haploid gametophyte pop
    
    dsex_starts = c(rep(1, asInteger(FtoM*dK)), rep(0, asInteger((1-FtoM)*dK)));
    p1.individuals.tag = dsex_starts;
"""

REPRO_ANGIO_DIO_P1 = """
    g_1 = genome1;
    g_2 = genome2;
    
    if (individual.tag == 1)
    {
        // determine how many ovules were fertilized, out of the total
        fertilizedOvules = rbinom(1, ovule_count, fertilization_rate);
        meiosis_reps = floor(fertilizedOvules/2);
        if (runif(1) <= Clone_rate)
            meiosis_reps = meiosis_reps*Clone_num;
        
        for (rep in 1:meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
        }
    
    }
    else //individual is male
    {
        meiosis_reps = floor(pollen_count/2);
        if (runif(1) <= Clone_rate)
            meiosis_reps = meiosis_reps*2;
        for (rep in 1:meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
        }
    }
"""

REPRO_ANGIO_DIO_P0  = """
    if (individual.tag == 1)  // females find male gametes to reproduce
    {
        if (pollen_comp == T)
        {
            pollen_pool = p0.sampleIndividuals(pollen_per_stigma, tag=0);   // sperm land on stigma
            for (pollen in pollen_pool)
            {
                pollen.setValue("fitness", p0.cachedFitness(pollen.index)); //store fitness value
                pollen.tag = 2;
            }
            
            if (pollen_pool.length()>0)
            {
            target_fitness = max(pollen_pool.getValue("fitness"));
            winners = pollen_pool[pollen_pool.getValue("fitness") == target_fitness];
            sperm = winners[0];
            }
            else sperm = p0.sampleIndividuals(1, tag=0);    // find a male
        }
        else
            sperm = p0.sampleIndividuals(1, tag=0); // find a male
        if (sperm.size() == 1)
        {
            child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
            sperm.tag = 2;
            
            if (runif(1) <= FtoM)
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
        fertilizedOvules = rbinom(1, ovule_count, fertilization_rate);
        meiosis_reps = floor(fertilizedOvules/2);
        if (runif(1) <= Clone_rate)
            meiosis_reps = meiosis_reps*Clone_num;
        
        for (rep in 1:meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 1;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 1;
        }
    
        meiosis_reps = floor(pollen_count/2);
        if (runif(1) <= Clone_rate)
            meiosis_reps = meiosis_reps*2;
        for (rep in 1:meiosis_reps)
        {
            breaks = sim.chromosome.drawBreakpoints(individual);
            p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL).tag = 0;
            p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL).tag = 0;
        }
"""

REPRO_PTER_HOMOSPORE_P0 = """
    if (runif(1) < gFtoM) // chance of making meristic (egg-bearing) gametophyte
    {
        reproduction_opportunity_count = 1;
        
        // clones give the focal individual extra opportunities to reproduce
        if (runif(1) <= gClone_rate)
            reproduction_opportunity_count = reproduction_opportunity_count + Clone_num;
        
        for (repro in seqLen(reproduction_opportunity_count))
        {
            if (runif(1) <= Self_rate)
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
                // this is selfing using two identical gametes – intragametophytic selfing
                // intergametophytic selfing might happen below, by chance

            else
            {
                sperm = p0.sampleIndividuals(1);
                
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                
                if (Maternal_weight > 0) //Mother's fitness affects sporophyte fitness; see survival()
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            
            }
        }
    }
"""

REPRO_PTER_HOMOSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;
    
    meiosis_reps = floor(Spore_num/2);
    reproduction_opportunity_count = 1;
    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= Clone_rate)
        meiosis_reps = meiosis_reps*Clone_num;
    
    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""

REPRO_PTER_HETEROSPORE_P0 = """
    if (runif(1) < gFtoM) // chance of making meristic (egg-bearing) gametophyte
    {
        reproduction_opportunity_count = 1;
        
        // clones give the focal individual extra opportunities to reproduce
        if (runif(1) <= gClone_rate)
            reproduction_opportunity_count = reproduction_opportunity_count + Clone_num;
        
        for (repro in seqLen(reproduction_opportunity_count))
        {
            if (runif(1) <= Self_rate)
                p1.addRecombinant(individual.genome1, NULL, NULL, individual.genome1, NULL, NULL);
                // this is selfing using two identical gametes – intragametophytic selfing
                // intergametophytic selfing might happen below, by chance

            else
            {
                sperm = p0.sampleIndividuals(1);
                
                child = p1.addRecombinant(individual.genome1, NULL, NULL, sperm.genome1, NULL, NULL);
                
                if (Maternal_weight > 0) //Mother's fitness affects sporophyte fitness; see survival()
                    child.setValue("maternal_fitness", subpop.cachedFitness(individual.index));
            
            }
        }
    }
"""

REPRO_PTER_HETEROSPORE_P1 = """
    g_1 = genome1;
    g_2 = genome2;
    
    meiosis_reps = floor(Spore_num/2);
    reproduction_opportunity_count = 1;
    // clones give the focal individual extra opportunities to reproduce
    if (runif(1) <= Clone_rate)
        meiosis_reps = meiosis_reps*Clone_num;
    
    for (rep in 1:meiosis_reps)
    {
        breaks = sim.chromosome.drawBreakpoints(individual);
        p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);
        p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);
    }
"""