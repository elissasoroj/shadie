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
    ACTIVATE, DEACTIVATE, EARLY_ANGIO, SURV,
    SUBSTITUTION, SUB_MUTS, EARLY1_ANGIO, REPRO_ANGIO_DIO_P1, 
    REPRO_ANGIO_DIO_P0, REPRO_ANGIO_MONO_P1,
    ANGIO_SURV_P0, LATE_ANGIO_DIO
)

DTYPES = ("d", "dio", "dioecy", "dioecious",)
MTYPES = ("m", "mono", "monoecy", "monecious",)

@dataclass
class SpermatophyteBase(ReproductionBase):
    lineage: str = field(default="Angiosperm", init=False)
    mode: str
    _file_in: str
    _chromosome: 'shadie.chromosome.ChromosomeBase'
    _sim_time: int
    _file_out: str

@dataclass
class Spermatophyte(SpermatophyteBase):
    """
    Reproduction mode based on angiosperms; appropriate for gymnosperms
    as well
    """
    spo_pop_size: int
    gam_pop_size: Union[None, int] = None
    microspore_pool: Union[None, int]=None,
    spo_mutation_rate: Union[None, float] = None
    gam_mutation_rate: float = 0.0
    spo_female_to_male_ratio: float.as_integer_ratio = (1,1)
    spo_clone_rate: float=0.0
    spo_clones_per:int = 1
    spo_flowers_per: int=5
    flower_ovules_per: int=15
    flower_pollen_per: int=100
    ovule_fertilization_rate: float=0.7
    pollen_success_rate: float=0.7
    pollen_comp: str="F"
    stigma_pollen_per: int=50
    spo_maternal_effect: float=0,
    spo_random_death_chance: float=0
    gam_random_death_chance: float=0

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        self.spo_female_to_male_ratio = (
            self.spo_female_to_male_ratio[0]/
            (self.spo_female_to_male_ratio[0]+self.spo_female_to_male_ratio[1]))

        if self.spo_mutation_rate:
            pass
        else:
            self.spo_mutation_rate = self.model.mutation_rate
            self.gam_mutation_rate = 0.0

        self.optimize()

        self.add_initialize_constants()
        self.add_early_haploid_diploid_subpops()   
        self.end_sim()     
        if self.mode in DTYPES:
            self.dioecious()
        elif self.mode in MTYPES:
            self.monoecious()
        else:
            raise ValueError(
                f"'mode' not recognized, must be in {DTYPES + MTYPES}")


    def optimize(self):
        """
        Optimizes sperm count to make simulation run faster; helps user
        check params
        """
        spo_microspores_generated = int(self.spo_pop_size
                    *(1-self.spo_female_to_male_ratio)
                    *self.spo_flowers_per
                    *self.flower_pollen_per
                    )
        spo_sperm_generated = spo_microspores_generated

        spo_megaspores_generated = int(self.spo_pop_size
                    *self.spo_female_to_male_ratio
                    *self.spo_flowers_per
                    *self.flower_ovules_per
                    )
        spo_eggs_generated = spo_megaspores_generated
        
        eggs_needed = self.spo_pop_size
        megaspores_needed = eggs_needed

        if self.pollen_comp == "T":
            sperm_needed = int(spo_eggs_generated*(self.stigma_pollen_per
                                    /self.flower_ovules_per))
            min_sperm_needed = int(eggs_needed*(self.stigma_pollen_per
                                        /self.flower_ovules_per))
        else:
            sperm_needed = spo_eggs_generated
            min_sperm_needed = eggs_needed

        optimized_megaspore_prop = eggs_needed/(eggs_needed+min_sperm_needed)

        if self.gam_pop_size is None:
            self.gam_pop_size = int(eggs_needed/optimized_megaspore_prop)
        else:
            suggested_gam_pop_size = int(eggs_needed/optimized_megaspore_prop)

        target_megaspore_prop_in_gam_pop = megaspores_needed/self.gam_pop_size

        megaspore_prop_check = spo_megaspores_generated/(
                    spo_megaspores_generated+spo_microspores_generated)

        logger.info("\nYour target megaspore:micospore proportion based "
            f"on `spo_pop_size` = {self.spo_pop_size} is: "
            f"{target_megaspore_prop_in_gam_pop}.\n"
            "\nYour calculatetd megaspore:micospore proportion is: "
            f"{megaspore_prop_check}.\n"
            f"This value is based on expectation of approximately " 
            f"{spo_eggs_generated} eggs and {spo_sperm_generated} "
            "sperm generated each generation."
            )

        eggs_alive_in_p0 =  int(megaspore_prop_check*self.gam_pop_size)
        
        if  eggs_alive_in_p0 < self.spo_pop_size:
            microspore_pool = sperm_needed
            total_spores = microspore_pool + spo_megaspores_generated
            new_megaspore_prop = spo_megaspores_generated/total_spores
            new_eggs_alive = int(new_megaspore_prop*self.gam_pop_size)
            self.microspore_pool = microspore_pool

            logger.info("\nGiven your current `gam_pop_size` = "
                f"{self.gam_pop_size}, your expected egg suvival is "
                f"{eggs_alive_in_p0}. This will not support `spo_pop_size`"
                f" of {self.spo_pop_size}.\n"
                f"Using `microspore_pool` = {self.microspore_pool} to "
                f"produce approximately "
                f"{self.microspore_pool} "
                f"sperm and {spo_eggs_generated} eggs per "
                f"generation, of which approximately {new_eggs_alive} "
                f"eggs will survive.\n"
                "New calculatetd megaspore:micospore proportion is "
                f"{new_megaspore_prop}\n"
                "If you are using the suggested gam_pop_size = "
                f"{suggested_gam_pop_size} "
                "and this new megaspore:micospore proportion is still not "
                "correct, then `flower_ovules_per` and `stigma_pollen_per` "
                "parameters likely need to be adjusted\n."
                )
        if self.gam_pop_size != suggested_gam_pop_size:
            logger.info("\nCurrent `gam_pop_size parameter = "
                f"{self.gam_pop_size}. Suggested min value is: "
                f"{suggested_gam_pop_size}.")

        if self.microspore_pool is None:
            pass

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
        constants["spo_clone_rate"] = self.spo_clone_rate
        constants["spo_clones_per"] = self.spo_clones_per
        constants["spo_flowers_per"] = self.spo_flowers_per
        constants["flower_ovules_per"] = self.flower_ovules_per
        constants["ovule_fertilization_rate"] = self.ovule_fertilization_rate
        constants["pollen_success_rate"] = self.pollen_success_rate
        constants["flower_pollen_per"] = self.flower_pollen_per
        constants["pollen_comp"] = self.pollen_comp
        constants["stigma_pollen_per"] = self.stigma_pollen_per
        constants["spo_maternal_effect"] = self.spo_maternal_effect
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
                scripts= EARLY1_ANGIO,
                comment="define Angiosperm subpops: diploid sporophytes, haploid gametophytes",
            )

    def end_sim(self):
        """
        adds late() call that ends the simulation and saves the .trees file
        """
        endtime = int(self._sim_time + 1)

        if self._file_in:
            ts_start = pyslim.load(self._file_in)
            sim_start = ts_start.max_root_time
            resched_end = int(endtime + sim_start)
            self.model.late(
                    time = resched_end, 
                    scripts = [
                    "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n",
                    f"sim.treeSeqOutput('{self._file_out}')"],
                    comment = "end of sim; save .trees file",
                )
        else:
            self.model.late(
                    time = endtime, 
                    scripts = [
                    "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n",
                    f"sim.treeSeqOutput('{self._file_out}')"],
                    comment = "end of sim; save .trees file",
                )


    def monoecious(self):
        """
        fills the script reproduction block with bryophyte-monoicous
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

        early_script = (
            EARLY_ANGIO.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'p0maternal_effect': "",
                'p1maternal_effect': "",
                'p0survival': ANGIO_SURV_P0}).lstrip())
        self.model.custom(survival_script)

        self.model.repro(
            population="p1",
            scripts=REPRO_ANGIO_MONO_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_ANGIO_DIO_P0,
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

    def dioecious(self):
        """
        fills the script reproduction block with bryophyte-dioicous
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

        early_script = (
            EARLY_ANGIO.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'p0maternal_effect': "",
                'p1maternal_effect': "",
                'p0survival': ""}).lstrip())
        self.model.custom(survival_script)

        self.model.repro(
            population="p1",
            scripts=REPRO_ANGIO_DIO_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_ANGIO_DIO_P0,
            comment="generates gametes from sporophytes"
            )

        substitution_str = ""
        for i in substitutions:
            substitution_str += "\n  ".join([i.strip(';') + ";\n    "])

        substitution_script = (
            SUBSTITUTION.format(**{'muts': substitution_str,
                'late': LATE_ANGIO_DIO}).lstrip())

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

        mod.reproduction.spermatophyte(
            mode='dio',
            spo_pop_size=1000, 
            gam_pop_size=4500,
            pollen_comp="T",
        )

    print(mod.script)
    #mod.run()
