#!/usr/bin/env python

"""
defaults shadiescripts
"""

from typing import Union

# 1 early()
SHADIE_POPS = """
	sim.addSubpop('p1', dK); //diploid sporophytes
	sim.addSubpop('p0', hK); //haploid gametophytes
}
"""

#early ()
EARLY_BRYO_DIO = """
{early()
{{
	if (sim.generation % 2 == 0) //diploids (p1) just generated haploid gametophytes
	{{
		//fitness affects gametophyte survival
		p0.fitnessScaling = (hK / p0.individualCount);
		
		// haploids get modified fitness, without dominance
		s1.active = 1;
		s2.active = 1;
		s3.active = 1;
		s4.active = 0;
		s5.active = 1;
		s6.active = 0;
	}}
	else //odd generations = gametophytes (p0) just generated sporophytes
	{{
		p1.fitnessScaling = dK / p1.individualCount; //fitness affects sporophytes
		// diploids get SLiM's standard fitness calculation, with dominance
		s1.active = 0;
		s2.active = 0;
		s3.active = 0;
		s4.active = 1;
		s5.active = 0;
		s6.active = 1;
	}}
}}

}
"""

SURV_BRYO ="""
{s3 survival(p1)
{{
	return F;
}}

s4 survival(p1)
{{
	//this code implements random death chance
	if (runif(1) < Death_chance)
		return F;
	else
		return NULL;
	
	// maternal effect
	maternal_effect = individual.getValue("maternal_fitness");
	
	if (!isNULL(maternal_effect))
	{{
		corrected_fitness = (maternal_effect * Maternal_weight) + fitness * (1 - Maternal_weight);
		return (draw < corrected_fitness);
	}}
	
	return NULL;
}}

s5 survival(p0) //even
{{
	//this code implements random death chance
	if (runif(1) < Death_chance)
		return F;
	else
		return NULL;
}}

s6 survival(p0) //odd
{{
	return F;
}}
}
"""

#no maternal effect in angio, BUT we could use same script above 
#with a default maternal_weight = 0 for angiosperm repro instead

SURV_ANGIO ="""
{s3 survival(p1)
{{
	return F;
}}

s4 survival(p1)
{{
	//this code implements random death chance
	if (runif(1) < Death_chance)
		return F;
	else
		return NULL;
}}

s5 survival(p0) //even
{{
	//this code implements random death chance
	if (runif(1) < Death_chance)
		return F;
	else
		return NULL;
}}

s6 survival(p0) //odd
{{
	return F;
}}
}
"""

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