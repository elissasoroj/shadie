#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""
from dataclasses import dataclass, field
from typing import Union
import pyslim
from loguru import logger
from shadie.reproduction.optimized.base_wf import ReproductionBase
from shadie.reproduction.optimized.scripts import (
    ACTIVATE, DEACTIVATE, EARLY, SURV,
    SUBSTITUTION, SUB_MUTS, 
    MATERNAL_EFFECT_P0, MATERNAL_EFFECT_P1,
    REPRO_PTER_HOMOSPORE_P1, REPRO_PTER_HOMOSPORE_P0,
    REPRO_PTER_HETEROSPORE_P1, REPRO_PTER_HETEROSPORE_P0,
    LATE_PTER_HOMOSPORE, LATE_PTER_HETEROSPORE
)

DTYPES = ("d", "dio", "dioicy", "dioicous", "heterosporous",)
MTYPES = ("m", "mono", "monoicy", "monoicous", "homosporous",)

@dataclass
class PteridophyteBase(ReproductionBase):
    lineage: str = field(default="Angiosperm", init=False)
    mode: str
    _file_in: str
    _chromosome: 'shadie.chromosome.ChromosomeBase'
    _sim_time: int
    _file_out: str

@dataclass
class Pteridophyte(PteridophyteBase):
    """
    Reproduction mode based on ferns and lycophytes
    """
    spo_pop_size: int
    gam_pop_size: int
    microspore_pool: Union[None, int]=None,
    spo_mutation_rate: Union[None, float] = None
    gam_mutation_rate: Union[None, float] = None
    spo_female_to_male_ratio: float.as_integer_ratio = (1,1)
    gam_female_to_male_ratio: float.as_integer_ratio = (1,1)
    spo_megaspores_per: int= 10
    spo_microspores_per: int=100
    gam_eggs_per_megaspore: int=1 
    gam_sperm_per_microspore: int=10
    spo_clone_rate: float=0.0
    spo_clones_per: float=1
    gam_clone_rate: float=0.0
    gam_clones_per: int=1
    spo_self_rate: float=0.0,
    gam_self_rate: float=0.0
    spo_maternal_effect: float=0
    gam_maternal_effect: float=0
    spo_random_death_chance: float=0
    gam_random_death_chance: float=0


    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        
        #calculate FtoM
        self.spo_female_to_male_ratio = (
            self.spo_female_to_male_ratio[0]/
            (self.spo_female_to_male_ratio[0]+self.spo_female_to_male_ratio[1]))

        #calculate gFtoM
        self.gam_female_to_male_ratio = (
            self.gam_female_to_male_ratio[0]/
            (self.gam_female_to_male_ratio[0]+self.gam_female_to_male_ratio[1]))

        #set up sporophyte and gametophyte mutation rates
        if self.spo_mutation_rate or self.gam_mutation_rate:
            assert self.spo_mutation_rate and self.gam_mutation_rate, (
                "You must define a mutation rate for both sporophyte "
                "and gametophyte generations.")
            if self.gam_mutation_rate:
                self.spo_mutation_rate = self.spo_mutation_rate
                self.gam_mutation_rate = self.gam_mutation_rate
        else:
            self.spo_mutation_rate = 0.5*self.model.mutation_rate
            self.gam_mutation_rate = 0.5*self.model.mutation_rate

        self.optimize()

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


    def optimize(self):
        spo_microspores_generated = int(self.spo_pop_size
                    *(1-self.spo_female_to_male_ratio)
                    *self.spo_microspores_per
                    )
        spo_sperm_generated = spo_microspores_generated*self.gam_sperm_per_microspore

        spo_megaspores_generated = int(self.spo_pop_size
                    *self.spo_female_to_male_ratio
                    *self.spo_megaspores_per
                    )
        spo_eggs_generated = spo_megaspores_generated*self.gam_eggs_per_megaspore
        eggs_needed = self.spo_pop_size
        megaspores_needed = eggs_needed/self.gam_eggs_per_megaspore

        target_megaspore_prop_in_gam_pop = megaspores_needed/self.gam_pop_size

        megaspore_prop_check = spo_megaspores_generated/(
                    spo_megaspores_generated+spo_microspores_generated)

        logger.info("\nYour target megaspore:micospore proportion based "
            f"on `gam_pop_size` = {self.gam_pop_size} is: "
            f"{target_megaspore_prop_in_gam_pop}.\n"
            f"If this value is very different from your supplied "
            f"`gam_female_to_male_ratio` = {self.gam_female_to_male_ratio}"
            f", then you may need to adjust these parameters."
            "\n\nYour calculatetd megaspore:micospore proportion is: "
            f"{megaspore_prop_check}.\n"
            f"This value is based on {spo_eggs_generated} eggs and "
            f"{spo_sperm_generated} sperm generated each generation"
            )

        eggs_alive_in_p0 =  int(megaspore_prop_check*self.gam_pop_size)
        
        if eggs_alive_in_p0 < self.spo_pop_size:
            microspore_pool = spo_megaspores_generated
            total_spores = microspore_pool + spo_megaspores_generated
            new_megaspore_prop = spo_megaspores_generated/total_spores
            new_eggs_alive = int(new_megaspore_prop
                        *self.gam_pop_size*self.gam_eggs_per_megaspore)
            self.microspore_pool = microspore_pool

            logger.info("\nGiven your current `gam_pop_size` = "
                f"{self.gam_pop_size}, your expected egg suvival is "
                f"{eggs_alive_in_p0}. This will not support `spo_pop_size`"
                f" of {self.spo_pop_size}.\n"
                f"Using `microspore_pool` = {self.microspore_pool} to "
                f"produce approximately "
                f"{self.microspore_pool*self.gam_sperm_per_microspore} "
                f"sperm and {spo_eggs_generated} eggs per "
                f"generation, of which approximately {new_eggs_alive} "
                f"eggs will survive.\n"
                "New calculatetd megaspore:micospore proportion is "
                f"{new_megaspore_prop}"
                ""
                )

    def add_initialize_constants(self):
        """
        Add defineConstant calls to init for new variables
        """
        constants = self.model.map["initialize"][0]['constants']
        constants["spo_pop_size"] = self.spo_pop_size
        constants["gam_pop_size"] = self.gam_pop_size
        constants["microspore_pool"] = self.microspore_pool
        constants["spo_mutation_rate"] = self.spo_mutation_rate
        constants["gam_mutation_rate"] = self.gam_mutation_rate
        constants["spo_female_to_male_ratio"] = self.spo_female_to_male_ratio
        constants["gam_female_to_male_ratio"] = self.gam_female_to_male_ratio
        constants["spo_megaspores_per"] = self.spo_megaspores_per
        constants["spo_microspores_per"] = self.spo_microspores_per
        constants["gam_eggs_per_megaspore"] = self.gam_eggs_per_megaspore
        constants["gam_sperm_per_microspore"] = self.gam_sperm_per_microspore
        constants["spo_clone_rate"] = self.spo_clone_rate
        constants["spo_clones_per"] = self.spo_clones_per
        constants["spo_self_rate"] = self.spo_self_rate
        constants["gam_self_rate"] = self.gam_self_rate
        constants["gam_clone_rate"] = self.gam_clone_rate
        constants["gam_clones_per"] = self.gam_clones_per
        constants["spo_maternal_effect"] = self.spo_maternal_effect
        constants["gam_maternal_effect"] = self.gam_maternal_effect
        constants["spo_random_death_chance"] = self.spo_random_death_chance
        constants["gam_random_death_chance"] = self.gam_random_death_chance


    def add_early_haploid_diploid_subpops(self):
        """
        add haploid and diploid life stages
        """
        if self._file_in:
            self.model.readfromfile()
        else:
            self.model.early(
                time=1,
                scripts= ["sim.addSubpop('p1', spo_pop_size)", 
                		"sim.addSubpop('p0', 0)",
                        "p1.individuals.tag=0"],
                comment="add p1, p0",
            )


    def end_sim(self):
        """
        adds late() call that ends the simulation and saves the .trees file
        """
        endtime = int(self._sim_time + 1)
        self.model.late(
                time = endtime, 
                scripts = [
                "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n",
                f"sim.treeSeqOutput('{self._file_out}')"],
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
        for mut in self._chromosome.mutations:
            i = i + 1
            idx = str("s" + str(i))
            active_script = ACTIVATE.format(**{'idx': idx}).lstrip()
            deactive_script = DEACTIVATE.format(**{'idx': idx}).lstrip()
            activate.append(active_script)
            deactivate.append(deactive_script)
            sub_MUTS = SUB_MUTS.format(**{'idx': idx, 'mut': mut}).lstrip()
            substitutions.append(sub_MUTS)
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
                'p1maternal_effect': MATERNAL_EFFECT_P1,
                'p0survival': ""}).lstrip())
        self.model.custom(survival_script)

        self.model.repro(
            population="p1",
            scripts=REPRO_PTER_HETEROSPORE_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_PTER_HETEROSPORE_P0,
            comment="generates gametes from sporophytes"
            )

        substitution_str = ""
        for i in substitutions:
            substitution_str += "\n  ".join([i.strip(';') + ";\n    "])

        substitution_script = (
            SUBSTITUTION.format(**{'muts': substitution_str,
                'late': LATE_PTER_HETEROSPORE}).lstrip())

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
        for mut in self._chromosome.mutations:
            i = i + 1
            idx = str("s" + str(i))
            active_script = ACTIVATE.format(**{'idx': idx}).lstrip()
            deactive_script = DEACTIVATE.format(**{'idx': idx}).lstrip()
            activate.append(active_script)
            deactivate.append(deactive_script)
            sub_muts = SUB_MUTS.format(**{'idx': idx, 'mut': mut}).lstrip()
            substitutions.append(sub_muts)
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
            SURV.format(**{'p0maternal_effect': MATERNAL_EFFECT_P0,
                'p1maternal_effect': MATERNAL_EFFECT_P1,
                'p0survival': ""}).lstrip())
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
            SUBSTITUTION.format(**{'muts': substitution_str,
                'late': LATE_PTER_HOMOSPORE}).lstrip())

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
            mode='homosporous',
            spo_pop_size=1000, 
            gam_pop_size=2000,
        )


    print(mod.script)
    #mod.run()
