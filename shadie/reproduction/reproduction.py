#!/usr/bin/env python

"""
Convenience functions for constructing each reproduction mode
"""

from typing import Union

#shadie defaults
SHADIE_POPS = """
1 early(){
	sim.addSubpop('p1', {diploid_pop}); //diploid sporophytes
	sim.addSubpop('p0', {haploid_pop}); //haploid gametophytes
}
"""


REPRO_BRYO = """
{reproduction(p1) {{
    g_1 = genome1;
	g_2 = genome2;
	//diploid undergoes meiosis exactly {meiosiscount}x
	for (meiosisCount in 1:{meiosiscount})
	{
		{malebreaks}

		{fembreaks}
}}

{reproduction(p0) {{  //creation of sporophyte from haploids
    if (runif(1) <= {cloning_rate})
    	clone_parent = p0.sampleIndividuals(1, tag = {female_tag});

    	child = p1.addCloned(clone_parent);
    	child.tag = {female_tag};

    if(runif(1) <= {haploid_selfing_rate})
    	self_parent = p0.sampleIndividuals(1, tag = {female_tag});

    	child = p1.addSelfed(self_parent);
    	child.tag = {female.tag};

    else
	    egg = p0.sampleIndividuals(1, tag = {female_tag});
	    sperm = p0.sampleIndividuals(1, tag = {male_tag});

	    child = p1.addRecombinant(egg.genome1, NULL, NULL, 
	        sperm.genome1, NULL, NULL);
	    //egg.tag = {used_tag};
	    //sperm.tag = {used_tag};
}}
"""

REPRO_MONO = """
//Hermaphroditic sporophytes
"""

REPRO_DIO = """
//Separate sex sporophytes
{reproduction(p1) {{
    g_1 = genome1;
	g_2 = genome2;
	//diploid undergoes meiosis exactly {meiosiscount}x
	for (meiosisCount in 1:{meiosiscount})
	{
		if (individual.getValue == {male_tag})
		{malebreaks}

		else if (individual.getValue == {female_tag})
		{fembreaks}
}}

{reproduction(p0) {{  //creation of sporophyte from haploids
    egg = p0.sampleIndividuals(1, tag = {female_tag});
    sperm = p0.sampleIndividuals(1, tag = {male_tag});

    child = p1.addRecombinant(egg.genome1, NULL, NULL, 
        sperm.genome1, NULL, NULL);
    if (runif(1) <= {haploid_ftom})
    	child.setValue = ("Sex", {female_tag});
	else
		child.tag = {male_tag};
    //egg.setValue = ("Sex", "{used_tag}");
}}


subpop.setValue("weights1", ifelse(has_m2, 2.0, 1.0));
"""

REPRO_HOMOSPORE = """
//Hermaphroditic gametophytes
"""

REPRO_HETEROSPORE = """
//Separate sex gametophytes

"""

SHADIE_FIT = """
fitness({mutation}) {{
    if (sim.generation % 2 == 0) //diploids (p1) just generated haploid gametophytes
        //gametophytes have no dominance effects
        return 1.0 + (mut.selectionCoeff * {haploid_fitness_scalar});
    else //odd generations = creation of diploids
        if (homozygous)
            return 1.0 + (mut.selectionCoeff * {diploid_fitness_scalar});
        else
            return 1.0 + (mut.mutationType.dominanceCoeff * mut.selectionCoeff
            * {diploid_fitness_scalar});
}}
"""

class Reproduction:

	def __init__(
		self, 
		mode = str(None),
		dipftom = float.as_integer_ratio(1/1),
		hapftom = float.as_integer_ratio(1/1),
		meiosiscount = int(5),
		lineage = string(None),
		dipfitscale = float(1.0),
		hapfitscale = float(1.0),
		dipselfrate = float(0.0),
		hapselfrate = float(0.0),
		):

		self.meiosiscount = meiosiscount
		self.dipftom = dipftom
		self.dipfitscale = dipfitscale
		self.gamfitscale = gamfitscale

		if meiosiscount/(dipftom[0] + dipftom[1]) is not int:
	    	logger.warning("female:male ratio must be compatible with meiosiscount")
	    else:
	    	pass

		self.totoffspring = meiosiscount/(dipftom[0] + dipftom[1])
	    self.females = dipftom[0]*totoffspring
	    self.males = dipftom[1]*totoffspring

        if mode = "d" or "dio" or "dioecy" or "dioecious" or "heterosporous":
        	self.dioecy()
        if mode = "m" or "mono" or "monoecy" or "monecious" or "homosporous":
        	self.monoecy()

	def dioecy(  
		self,     
	    femtag  = "female",
	    maletag = int(2),
	    hermtag  = int(3),
	    usedtag = int(5),
	    ):
	    """
	    Sets up a dioecious reproduction mode
	    """

	    #write the script
	    self.script[('reproduction', None)] = (
	            REPRO_DIO.format(**{
	                "female_tag": femtag,
	                "male_tag": maletag,
	                "hermaphrodite_tag": hermtag,
	                "used_tag": usedtag,
	                "scripts": "\n  ".join([i.strip(";") + ";" for i in scripts]),
	            })
	        )

	def breaks():
	"""
	Breaks for p1 (diploid) reproduction callback:
	"""
	    fembreaks =  ""
	    if females/2 is int:
	    	for i in range(1, females):
	    		fembreaks.append("breaks = sim.chromosome.drawBreakpoints(individual);"
	    			f"f_{i} = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);"
	    			f"f_{i}.tag = {femtag};"
	    			)
	    else:
	    	for i in range(1, females-1):
	    		fembreaks.append("breaks = sim.chromosome.drawBreakpoints(individual);"
	    			f"f_{i} = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);"
	    			f"f_{i}.tag = {femtag};"
	    			)
	    	#last block for odd  number of female offspring
	    	fembreaks.append("breaks = sim.chromosome.drawBreakpoints(individual);"
	    		"if (runif(1) <= 0.5)"
					f"f_{females} = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);"
				"else"
					f"f_{females} = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);"
				f"f_{females}.tag = {femtag};"
			)

		malebreaks =  ""
	    if males/2 is int:
	    	for i in range(1, males):
	    		malebreaks.append("breaks = sim.chromosome.drawBreakpoints(individual);"
	    			f"f_{i} = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);"
	    			f"f_{i}.tag = {maletag};"
	    			)
	    else:
	    	for i in range(1, males-1):
	    		malebreaks.append("breaks = sim.chromosome.drawBreakpoints(individual);"
	    			f"f_{i} = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);"
	    			f"f_{i}.tag = {maletag};"
	    			)
	    	#last block for odd  number of female offspring
	    	malebreaks.append("breaks = sim.chromosome.drawBreakpoints(individual);"
	    		"if (runif(1) <= 0.5)"
					f"f_{males} = p0.addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL);"
				"else"
					f"f_{males} = p0.addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL);"
				f"f_{males}.tag = {maletag};"
			)


	def monoecy(
		hermtag =  int(3),
		usedtag = int(5),
		meiosiscount = int()
		):5
	    """
	    Sets up a monoecious reproduction  mode
	    """
	    self.femtag = hermtag
	    self.maletag = hermtag


	def bryophytes(
		hap_pop = int(1000),
		dip_pop = int(1000),
		):
		"""
	    Reproduction mode based on mosses
	    """
	    #early callback:

	    rpdn_early = """
	    if (sim.generation % 2 == 0) //diploids (p1) just generated gametophytes
		{
			//fitness affects gametophyte survival:
			p0.fitnessScaling = {hap_pop}/p0.individualCount;
			p1.fitnessScaling = 0.0; // diploids die after producing spores
			for (ind in p0.individuals)
				if (ind.tag == {used_tag})
				ind.tag == {female_tag};
		}
		
		else //odd generations = gametophytes (p0) just generated sporophytes
		{	
			p0.fitnessScaling = 0.0; //gametophytes die after producing sporophytes
			//fitness affects sporophyte survival:
			p1.fitnessScaling = {dip_pop}/ p1.individualCount; 
		}
	    """

	    #at some point this needs to be appended to scripts:
	    rpdndict = {(early, None), rpdn_early}

	def lycophytes():
		"""
	    Reproduction mode based on lycophytes
	    """

	def monilophytes():
	 	"""
	    Reproduction mode based on ferns
	    """

	def angiosperm():
		"""
	    Reproduction mode based on flowering  plants
	    """



