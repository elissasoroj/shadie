#!/usr/bin/env python

"""
Produces reproduction scripts
"""

#package imports

class Lifecycle:
    "creates reproduction callbacks for SLiM script"

	def __init__(self):


	def WF(self):
		"standard WF model, implemented in base SLiM (diploid hermaphrodite)"
		rpdn = (
			"reproduction() {\n"
			f"K = sim.getValue('{self.Ne}'');\n"
			"for (i in seqLen(K))\n{"
			f"firstParent = {self.rpdndict[self.demog.loc[0]['src']]}.sampleIndividuals(1);\n"
			f"secondParent = {self.rpdndict[self.demog.loc[0]['src']]}.sampleIndividuals(1);\n"
			f"{self.rpdndict[self.demog.loc[0]['src']]}.addCrossed(firstParent, secondParent);\n"
			"}\nself.active = 0;"
			"}")

	def bryophyte(self):
		pass

	def pteridophyte(self):
		""
		rpdn = (
            "\nreproduction() {\n"
            "g_1 = genome1;\n"
            "g_2 = genome2;\n"
            "for (meiosisCount in 1:5)\n"
            "{\n"
                "if (individual.sex == 'M')\n"
                "{\n"
                    "breaks = sim.chromosome.drawBreakpoints(individual);\n"
                    f"s_1 = {self.rpdndict[self.demog.loc[0]['src']]}."
                    "addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL, 'M');\n"
                    f"s_2 = {self.rpdndict[self.demog.loc[0]['src']]}."
                    "addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL, 'M');\n"
                    "\n"
                    "breaks = sim.chromosome.drawBreakpoints(individual);\n"
                    f"s_3 = {self.rpdndict[self.demog.loc[0]['src']]}."
                    "addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL, 'M');\n"
                    f"s_4 = {self.rpdndict[self.demog.loc[0]['src']]}."
                    "addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL, 'M');\n"
                "}\n"
                "else if (individual.sex == 'F')\n"
                "{\n"
                    "breaks = sim.chromosome.drawBreakpoints(individual);\n"
                    "if (runif(1) <= 0.5)\n"
                        f"e = {self.rpdndict[self.demog.loc[0]['src']]}."
                        "addRecombinant(g_1, g_2, breaks, NULL, NULL, NULL, 'F');\n"
                    "else\n"
                        f"e = {self.rpdndict[self.demog.loc[0]['src']]}."
                        "addRecombinant(g_2, g_1, breaks, NULL, NULL, NULL, 'F');\n"
                "}\n"
            "}\n"
            "}\n"

            f"reproduction({self.rpdndict[self.demog.loc[0]['src']]}, 'F')\n"
            "{\n"
                f"mate = {self.rpdndict[self.demog.loc[0]['src']]}.sampleIndividuals(1, sex='M', tag=0);\n"
                " mate.tag = 1;"
                
                f"child = {self.demog.loc[0]['src']}.addRecombinant(individual.genome1, "
                "NULL, NULL, mate.genome1, NULL, NULL);\n"
            "}\n"
            "early()\n"
            "{\n"
                "if (sim.generation % 2 == 0)\n"
                "{\n"
                    f"{self.demog.loc[0]['src']}.fitnessScaling = 0.0;\n"
                    f"{self.rpdndict[self.demog.loc[0]['src']]}.individuals.tag = 0;\n"
                    "sim.chromosome.setHotspotMap(0.0);\n"
                "}\n"
                "else\n"
                "{\n"
                    f"{self.rpdndict[self.demog.loc[0]['src']]}.fitnessScaling = 0.0;\n"
                    f"{self.demog.loc[0]['src']}.fitnessScaling = {self.Ne} / {self.demog.loc[0]['src']}.individualCount;\n"
                    f"sim.chromosome.setHotspotMap(1.0);\n"
                "}\n"
            "}\n"
            )

        fit1 = "" 
        for key in self.mutationlist.mutationdict:
            fitness1 =+ (
                f"fitness({key})" + "{\n"
                "if (sim.generation % 2 == 0)\n" #creation of haploids occurs in even gens
                    "return 1.0 + mut.selectionCoeff;\n"
                "else\n" #odd generations = creation of diploids happens in odd gens
                    "if (homozygous)\n"
                        "return 1.0 + mut.selectionCoeff;\n"
                    "else\n"
                        "return 1.0 + mut.mutationType.dominanceCoeff * mut.selectionCoeff;\n"
                "}"
            )


    def gymnosperm(self):
		pass

	def angiosperm(self):
		pass


    def clonalhap(self):

        rpdn1 = (
            "reproduction() {\n"
            f"{self.rpdndict[self.demog.loc[0]['src']]}.addRecombinant("
            "genome1, NULL, NULL, NULL, NULL, NULL);\n"
            "\n"

            "early() {\n"
            f"{self.rpdndict[self.demog.loc[0]['src']]}.fitnessScaling = \n"
            f"{self.Ne} / {self.rpdndict[self.demog.loc[0]['src']]}.individualCount;\n"
            "\n"
           )

        fit1 = "" 
        for key in self.mutationlist.mutationdict:
            fitness1 =+ (f"fitness({key}) " + "{return 1.0 + mut.selectionCoeff;}")

