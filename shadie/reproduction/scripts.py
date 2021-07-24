#!/usr/bin/env python

"""
defaults shadiescripts
"""

# 1 early()
SHADIE_POPS = """
	sim.addSubpop('p1', dK); //diploid sporophytes
	sim.addSubpop('p0', hK); //haploid gametophytes
}
"""

#early ()
EARLY = """
{early()
{{
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

}
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
// maternal effect
	maternal_effect = individual.getValue("maternal_fitness");
	
	if (!isNULL(maternal_effect))
	{{
		corrected_fitness = (maternal_effect * Maternal_weight) + fitness * (1 - Maternal_weight);
		return (draw < corrected_fitness);
	}}
	
	return NULL;
"""

DEATH = """
return F;
"""

#------
#SHADIE SURVIVAL BLOCKS

SURV ="""
{s1 survival(p1)
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
}
"""

#-----------------------------------------------
#late() **for every mut!
LATE = """
{{
if (sim.generation % 2 == 0) //gametophytes have just undergone fitness selection
{{
	mut{idx} = sim.mutationsOfType({mut});
	freq{idx} = sim.mutationFrequencies(NULL, mut_{mut});
	if (any(freq{idx} == 0.5))
		sim.subpopulations.genomes.removeMutations(mut_{mut}[freq{idx} == 0.5], T);
}}
}}
"""