#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""
from typing import Union
from dataclasses import dataclass, field
from shadie.reproduction.base import ReproductionBase
from shadie.reproduction.base_scripts import (
    ACTIVATE, DEACTIVATE, EARLY, SURV,
    SUBSTITUTION, SUB_INNER, MATERNAL_EFFECT,
    REPRO_PTER_HOMOSPORE_P1, REPRO_PTER_HOMOSPORE_P0,
    REPRO_PTER_HETEROSPORE_P1, REPRO_PTER_HETEROSPORE_P1
)

DTYPES = ("d", "dio", "dioicy", "heterosporous", "dioicous")
MTYPES = ("m", "mono", "monoicy", "homosporous", "monoicous")

@dataclass
class PteridophyteBase(ReproductionBase):
    lineage: str = field(default="Angiosperm", init=False)
    mode: str
    chromosome: 'shadie.chromosome.ChromosomeBase'
    simtime: int
    fileout: str

@dataclass
class Pteridophyte(PteridophyteBase):
    """
    Reproduction mode based on ferns and lycophytes
    """
    diploid_ne: int
    haploid_ne: int
    female_to_male_ratio: float.as_integer_ratio = (1,1)
    gam_female_to_male_ratio: float.as_integer_ratio = (1,1)
    clone_rate: float=0.0
    gam_clone_rate: float=0.0
    spores_per_sporophyte: int=100
    selfing_rate: float=0.0
    maternal_effect_weight: float=0
    random_death_chance: float=0

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        self.female_to_male_ratio = (
            self.female_to_male_ratio[1]/
            (self.female_to_male_ratio[0]+self.female_to_male_ratio[1]))

        self.gam_female_to_male_ratio = (
            self.gam_female_to_male_ratio[1]/
            (self.gam_female_to_male_ratio[0]+self.gam_female_to_male_ratio[1]))

        self.add_initialize_constants()
        self.add_early_haploid_diploid_subpops() 
        self.end_sim()       
        if self.mode in DTYPES:
            self.heterosporous()
        elif self.mode in MTYPES:
            self.homosporous()
        else:
            raise ValueError(
                f"'mode' not recognized, must be in {DTYPES + MTYPES}")


    def add_initialize_constants(self):
        """
        Add defineConstant calls to init for new variables
        """
        constants = self.model.map["initialize"][0]['constants']
        constants["dK"] = self.diploid_ne
        constants["hK"] = self.haploid_ne
        constants["Death_chance"] = self.random_death_chance
        constants["FtoM"] = self.female_to_male_ratio
        constants["gFtoM"] = self.gam_female_to_male_ratio
        constants["Clone_rate"] = self.clone_rate
        constants["gClone_rate"] = self.gam_clone_rate
        constants["Spore_num"] = self.spores_per_sporophyte
        # constants["Clone_num"] = self.clone_number
        constants["Maternal_weight"] = self.maternal_effect_weight
        constants["Self_rate"] = self.selfing_rate


    def add_early_haploid_diploid_subpops(self):
        """
        add haploid and diploid life stages
        """
        self.model.early(
            time=1,
            scripts= ["sim.addSubpop('p1', dK)", "sim.addSubpop('p0', hK)"],
            comment="add p1, p0",
        )


    def end_sim(self):
        """
        adds late() call that ends the simulation and saves the .trees file
        """
        endtime = int(self.simtime + 1)

        self.model.late(
                time = endtime, 
                scripts = [
                #"sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n",
                f"sim.treeSeqOutput('{self.fileout}')"],
                comment = "end of sim; save .trees file",
            )

    def heterosporous(self):
        """
        Model is not ready yet
        """
        # fitness callback:
        i = 4
        activate = []
        deactivate = []
        substitutions = []
        for mut in self.chromosome.mutations:
            i = i + 1
            idx = str("s" + str(i))
            active_script = ACTIVATE.format(**{'idx': idx}).lstrip()
            deactive_script = DEACTIVATE.format(**{'idx': idx}).lstrip()
            activate.append(active_script)
            deactivate.append(deactive_script)
            sub_inner = SUB_INNER.format(**{'idx': idx, 'mut': mut}).lstrip()
            substitutions.append(sub_inner)
            self.model.fitness(
                idx=idx,
                mutation=mut,
                scripts="return 1 + mut.selectionCoeff",
                comment="gametophytes have no dominance effects",
            )
            
        activate_str = ""
        deactivate_str = ""
        for i in activate:
            activate_str += "\n  ".join([i.strip(';') + ";\n    "])

        for i in deactivate:
            deactivate_str += "\n  ".join([i.strip(';') + ";\n    "])

        self.active = activate_str

        early_script = (
            EARLY.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'p0maternal_effect': MATERNAL_EFFECT,
                'p1maternal_effect': "",
                'p0survival': "return NULL;"}).lstrip())
        self.model.custom(survival_script)

        self.model.repro(
            population="p1",
            scripts=REPRO_PTER_HOMOSPORE_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_PTER_HOMOSPORE_P0,
            comment="generates gametes from sporophytes"
            )

        substitution_str = ""
        for i in substitutions:
            substitution_str += "\n  ".join([i.strip(';') + ";\n    "])

        substitution_script = (
            SUBSTITUTION.format(**{'inner': substitution_str}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
            )

    def homosporous(self):
        """
        fills the script reproduction block with pteridophyte-homosporous
        """

        # fitness callback:
        i = 4
        activate = []
        deactivate = []
        substitutions = []
        for mut in self.chromosome.mutations:
            i = i + 1
            idx = str("s" + str(i))
            active_script = ACTIVATE.format(**{'idx': idx}).lstrip()
            deactive_script = DEACTIVATE.format(**{'idx': idx}).lstrip()
            activate.append(active_script)
            deactivate.append(deactive_script)
            sub_inner = SUB_INNER.format(**{'idx': idx, 'mut': mut}).lstrip()
            substitutions.append(sub_inner)
            self.model.fitness(
                idx=idx,
                mutation=mut,
                scripts="return 1 + mut.selectionCoeff",
                comment="gametophytes have no dominance effects",
            )
            
        activate_str = ""
        deactivate_str = ""
        for i in activate:
            activate_str += "\n  ".join([i.strip(';') + ";\n    "])

        for i in deactivate:
            deactivate_str += "\n  ".join([i.strip(';') + ";\n    "])

        self.active = activate_str

        early_script = (
            EARLY.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'p0maternal_effect': "",
                'p1maternal_effect': MATERNAL_EFFECT,
                'p0survival': "return NULL;"}).lstrip())
        self.model.custom(survival_script)

        self.model.repro(
            population="p1",
            scripts=REPRO_PTER_HOMOSPORE_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_PTER_HOMOSPORE_P0,
            comment="generates gametes from sporophytes"
            )

        substitution_str = ""
        for i in substitutions:
            substitution_str += "\n  ".join([i.strip(';') + ";\n    "])

        substitution_script = (
            SUBSTITUTION.format(**{'inner': substitution_str}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
            )


if __name__ == "__main__":


    import shadie
    with shadie.Model() as mod:
        
        # define mutation types
        m0 = shadie.mtype(0.5, 'n', 0, 0.4)
        m1 = shadie.mtype(0.5, 'g', 0.8, 0.75)
        #I suggest we add a checkpoint that calculates the average
        #fitness of mutations input by the user. If fitness is too high
        #the simuulation will lag tremendously. 
        # OK: a good use case for logger.warning('fitness is too high...')
        
        # define elements types
        e0 = shadie.etype([m0, m1], [1, 2])
        e1 = shadie.etype([m1], [1])
        
        # design chromosome of elements
        chrom = shadie.chromosome.random(
            genome_size=20000,
            noncds=e0,
            intron=e0,
            exon=e1,
        )

        # init the model
        mod.initialize(chromosome=chrom)

        mod.reproduction.pteridophyte(
            mode='mono',
            chromosome = chrom,
            diploid_ne=1000, 
            haploid_ne=1000,
        )


    print(mod.script)
    #mod.run()
