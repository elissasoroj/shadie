#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import subprocess
import numpy as np
import pandas as pd
import pyslim
from loguru import logger

from shadie.chromosome import Chromosome
from shadie.demography import Demography

class Shadie(object):
    "Produces a shadie.slim script for running a simulation in SLiM3"

    def __init__(
        self,
        tree=None,              # reads in 
        chromosome = None,      # reads in chromosome object from Chromosome class
        Ne = 2000,             # K = carrying capacity
        nsamples=2,             # number of sampled haplotypes per tip in final data 
        reproduction = None,    # defines how gametes get selected and replicate
        recomb = 1e-9,          # sets rate in `initializeRecombinationRate`, also accepts map
        generations = 10000,     #required if no tree object is supplied
        model = "WF"
        ):
        """
        Builds script to run SLiM3 simulation

        Parameters:
        -----------
        tree: (str)
            Optional. A newick string or Toytree object of a species tree with edges in 
            SLiM 'generation' units

        Ne (int): default = 1000
            The effective population size. This value will be set to all edges of the tree.

        organism (str):
            Options: "pteridophyte", "bryophyte", "angiosperm".
            Defines the haploid/diploid lifecycle and how selection will act at different
            life stages. Also defines how individuals replicate and hoow gametes are generated

        recomb (float):
            The per-site per-generation recombination rate.

        """
        self.model = model         # nonWF is needed for repoduction 
        self.recomb = recomb
        self.Ne = Ne
        self.chromosome = chromosome
        self.tree = tree
        self.reproduction = reproduction
        self.gens = generations
        
        if self.tree is not None:
            initdemog = Demography(self.tree)
            initdemog.get_demog()
            self.demog = initdemog.demog.sort_values("gen")
            self.treeheight = int(tree.treenode.height)
            self.rpdndict = initdemog.rpdndict
        else:
            defdemog = [{'gen': "1", 'src': 'p1', 'Ne': self.Ne}]
            self.demog = pd.DataFrame(defdemog)
            self.rpdndict = {'p1': 'p0'}
            self.treeheight = self.gens
            logger.warning("if no tree is provided, 'generations' "
                "argument must be provided (defines length of simulation in "
                f"generations. Current value = {self.gens}")


        if isinstance(self.chromosome, Chromosome): 
            self.mutationlist = self.chromosome.mutationlist
            self.elementlist = self.chromosome.elementlist
            self.genome = self.chromosome.genome
            self.gensize = self.chromosome.gensize

        elif self.chromosome is None:
            gene = Chromosome()
            self.genome = gene.genome 
            self.mutationlist = gene.mutationlist
            self.elementlist = gene.elementlist
            self.genome = gene.genome
            self.gensize = gene.gensize

        else:
            raise ValueError("please input valid Chromosome class object")

        exonstart = []
        exonstop = []

        
        for index, row in self.genome.iterrows():
            if row["type"] == "exon":
                exonstart.append(row['start'])
                exonstop.append(row['finish'])

        genemap = pd.DataFrame(list(zip(exonstart, exonstop)),
              columns=['exonstart', 'exonstop'])

        self.genemap = genemap


    def write(self, filename="shadie.slim"):
        "writes the .slim script; optional to provide filename as 'filename.slim'"
        
        self.filename = filename #add check .slim extension on filename
        outname = filename[:-5] + ".trees"
        self.outname = outname
        #else:     

        #write initialize callbacks
        script = open(self.filename, "w") #overwrites old file
        init1 = (
            "initialize() {\ninitializeSLiMModelType('"+self.model+"');\n"
            "initializeSLiMOptions(nucleotideBased=T);\n"
            "initializeHotspotMap(1.0);\n"
            "initializeSex('A');\n"
            "initializeTreeSeq();\n"
            f"initializeAncestralNucleotides(randomNucleotides({int(self.gensize)}));\n")
        init2 = ""
        for key in self.mutationlist.mutationdict:
            init2 += f"initializeMutationTypeNuc({self.mutationlist.mutationdict[key]});\n"
            f"{key}.convertToSubstitution = T;\n"
        
        init3 = ""
        for key in self.elementlist.elementdict:
            init3 += f"initializeGenomicElementType({self.elementlist.elementdict[key]});\n"

        init4 = ""
        for index, row in self.genome.iterrows():
            init4 += f"initializeGenomicElement({row['eltype']}, {int(row['start'])}, {int(row['finish'])});\n"

        initfinal = "initializeRecombinationRate("+str(self.recomb)+");\n}\n"

        #######

        #write the first line:
        if self.reproduction is not None:   
            rpdn0 = (f"sim.addSubpop('{self.rpdndict[self.demog.loc[0]['src']]}', 0);\n")
        else:
            rpdn0 =  ""

        start = (
            "\n1 early(){\n"
            f"sim.addSubpop('{self.demog.loc[0]['src']}', {self.demog.loc[0]['Ne']});\n"
            + rpdn0 + 
            "}\n"
            )
        #write the reproduction callbacks
        if self.reproduction is not None:
            rep1 = (
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
        else:
            rep1 = ""

        #write the demography
        if self.tree is not None:
            gens = ""
            for index, row in self.demog.iterrows():
                gens += (
                    f"{row['gen']}" + "{\n"
                    f"sim.addSubpopSplit('{row['child0']}', {row['Ne']}, {row['src']});\n"
                    f"sim.addSubpopSplit('{row['child1']}',{row['Ne']}, {row['src']});\n"
                    f"{row['src']}.setSubpopulationSize(0);" + "\n}\n"
                    )
        else:
            gens = ""

        final = (
            f"{self.treeheight} late()" + "{\n"
            f"sim.treeSeqOutput('{self.outname}');\n"
            f"sim.outputFull('model_output.txt');\n"
            "}"
            )

        script.write(init1 + init2 + init3 + init4 + initfinal + start + gens + rep1 + final)
        script.close

    def run(self):
        "calls SLiM to run the simulations"
        # Run the SLiM model and load the resulting .trees
        
        subprocess.check_output(["slim", "-m", "-s", "0", self.filename])
        
    def postsim(self):
        "post-SLiMulation analysis"
        self.ts = pyslim.load(self.outname)

        M = [[0 for _ in pyslim.NUCLEOTIDES] for _ in pyslim.NUCLEOTIDES]
        for mut in self.ts.mutations():
            mut_list = mut.metadata["mutation_list"]
            k = np.argmax([u["slim_time"] for u in mut_list])
            derived_nuc = mut_list[k]["nucleotide"]
            if mut.parent == -1:
                acgt = self.ts.reference_sequence[int(mut.position)]
                parent_nuc = pyslim.NUCLEOTIDES.index(acgt)
            else:
                parent_mut = self.ts.mutation(mut.parent)
                assert parent_mut.site == mut.site
                parent_nuc = parent_mut.metadata["mutation_list"][0]["nucleotide"]
            M[parent_nuc][derived_nuc] += 1
        
        counts = 0        
        print("{}\t{}\t{}".format('ancestr', 'derived', 'count'))
        for j, a in enumerate(pyslim.NUCLEOTIDES):
            for k, b in enumerate(pyslim.NUCLEOTIDES):
                counts += M[j][k]
                print("{}\t{}\t{}".format(a, b, M[j][k]))

        print(
            f"\nNumber of mutations: {counts}\n"
            " ---------------------\n"
            "Simulation settings\n\n"
            f"Carrying Capacity: {self.Ne}\n"
            f"Generations: {self.gens}\n"
            f"Tree: {self.tree}\n"
            f"Reproduction: {self.reproduction}\n"
            )


    def reproduction(self):
        "defines haploid/diploid life cycle offspring generation, selection, fitness effects"

        if self.reproduction == "pter":
            #define reproduction - sporophyte produces male or hermaphroditic gametophytes
            #define sex ratio?  - this may be relevant (e.g. hermaphroodite gametophytes 
            #induce formation of male-only gametophyte around them)
            #define fitness effects - selection on gametophyte
            #**no selfing with male/female??
            pass

        elif self.reproduction == "bryo":
            #define reproduction - male (sporophyte only?), female
            #define sex ratio?  - females live longer than males, 
            #define fitness effects - selection on gametophyte occurrs for longer than selection 
                # on sporophyte
            #**no selfing with male/female??
            pass

        elif self.reproduction == "angio":
            #define reproduction - hermaphrodite spermatophytes, male/female gameotphyte
            #define sex ratio?  - can sex ratio be set for haploid? male >> female for gametophyte
            #define fitness effects - selection mainly acts on sporophyte
            pass

        else:
            #base SLiM hermaphrodite? 
            pass



if __name__ == "__main__":
    import toytree

    #Make the tree
    tree = toytree.rtree.unittree(ntips=10, treeheight=1e4, seed=123)
    randtree = tree.set_node_values(
        feature="Ne", 
        values={i: np.random.randint(10000, 100000) for i in tree.idx_dict}
    )
    tree_sim = Shadie(tree = randtree)
