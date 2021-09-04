#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""
import pyslim
from loguru import logger
from typing import Union
from dataclasses import dataclass, field
from shadie.reproduction.optimized.base_wf import ReproductionBase
from shadie.reproduction.optimized.scripts import (
    ACTIVATE, DEACTIVATE, EARLY, SURV, MATERNAL_EFFECT_P0,
    MATERNAL_EFFECT_P1, SUBSTITUTION, SUB_MUTS, REPRO_BRYO_DIO_P1, 
    REPRO_BRYO_DIO_P0, REPRO_BRYO_MONO_P1, REPRO_BRYO_MONO_P0,
    LATE_BRYO_DIO, LATE_BRYO_MONO
)

DTYPES = ("d", "dio", "dioicy", "dioicous", "heterosporous",)
MTYPES = ("m", "mono", "monoicy", "monoicous", "homosporous",)

@dataclass
class BryophyteBase(ReproductionBase):
    lineage: str = field(default="Bryophyte", init=False)
    mode: str
    _file_in: str
    _chromosome: 'shadie.chromosome.ChromosomeBase'
    _sim_time: int
    _file_out: str

@dataclass
class Bryophyte(BryophyteBase):
    """
    Reproduction mode based on mosses, hornworts, and liverworts
    """
    spo_pop_size: int
    gam_pop_size: int
    spo_mutation_rate: Union[None, float] = None
    gam_mutation_rate: Union[None, float] = None
    gam_female_to_male_ratio: float.as_integer_ratio = (1,1)
    spo_megaspores_per: int= 10
    spo_microspores_per: int=100
    gam_eggs_per_megaspore: int=1 
    gam_sperm_per_microspore: int=10
    gam_clone_rate: float=0.0
    gam_clone_number: int=1
    spo_self_rate: float=0
    gam_self_rate: float=0
    gam_maternal_effect: float=0
    spo_random_death_chance: float=0
    gam_random_death_chance: float=0


    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        
        #calculate gFtoM
        self.gam_female_to_male_ratio = (
            self.gam_female_to_male_ratio[0]/
            (self.gam_female_to_male_ratio[0]+self.gam_female_to_male_ratio[1])
            )

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

        #optimize spore numbers
        eggs_per_gen = (self.gam_pop_size*self.gam_female_to_male_ratio*
                    (1-self.gam_random_death_chance)*
                    self.spo_megaspores_per*self.gam_eggs_per_megaspore)
        
        sperm_per_gen = (self.gam_pop_size*
            (1-self.gam_female_to_male_ratio)*self.spo_microspores_per*
            (1-self.gam_random_death_chance)*self.gam_sperm_per_microspore)

        fertilization_chance = sperm_per_gen/eggs_per_gen

        eggs_alive = self.gam_pop_size*(eggs_per_gen/sperm_per_gen)
        
        logger.info("With these simulation parameters, you will "
            f"generate {int(eggs_per_gen)} eggs and {int(sperm_per_gen)} "
            "sperm each generation. Likelihood of fertilization is "
            f"{min(int(fertilization_chance*100), 100)}%. Use `.optimize()` to choose "
            "optimized parameter values.\n"
            "You are generating enough eggs to reach spo_pop_size: "
            f"{eggs_per_gen > self.spo_pop_size}\n"
            "Enough eggs will live to next gen: "
            f"{eggs_alive > self.spo_pop_size}\n"
            )
        if not eggs_alive > self.spo_pop_size:
            logger.info("your egg:sperm ratio is too low given "
                "gam_pop_size. \nYour expected egg suvival is "
                f"{int(eggs_alive)} eggs per gen"
                )

        self.eggs_per_gen = eggs_per_gen
        self.sperm_per_gen = sperm_per_gen
        self.fertilization_chance = fertilization_chance

        self.add_initialize_constants()
        self.add_early_haploid_diploid_subpops() 
        self.end_sim()        
        if self.mode in DTYPES:
            self.dioicous()
        elif self.mode in MTYPES:
            self.monoicous()
        else:
            raise ValueError(
                f"'mode' not recognized, must be in {DTYPES + MTYPES}")

    def add_initialize_constants(self):
        """
        Add defineConstant calls to init for new variables
        """
        constants = self.model.map["initialize"][0]['constants']
        constants["spo_pop_size"] = self.spo_pop_size
        constants["gam_pop_size"] = self.gam_pop_size
        constants["spo_mutation_rate"] = self.spo_mutation_rate
        constants["gam_mutation_rate"] = self.gam_mutation_rate
        constants["gam_female_to_male_ratio"] = self.gam_female_to_male_ratio
        constants["spo_megaspores_per"] = self.spo_megaspores_per
        constants["spo_microspores_per"] = self.spo_microspores_per
        constants["gam_eggs_per_megaspore"] = self.gam_eggs_per_megaspore
        constants["gam_sperm_per_microspore"] = self.gam_sperm_per_microspore
        constants["gam_clone_rate"] = self.gam_clone_rate
        constants["gam_clone_number"] = self.gam_clone_number
        constants["spo_self_rate"] = self.spo_self_rate
        constants["gam_self_rate"] = self.gam_self_rate
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
                scripts=["sim.addSubpop('p1', spo_pop_size)",
                        "sim.addSubpop('p0', 0)",
                        "p1.individuals.tag=0"],
                comment="define Bryophyte subpops: diploid sporophytes, haploid gametophytes",
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

    def dioicous(self):
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
            EARLY.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'p0maternal_effect': MATERNAL_EFFECT_P0, #affects p1
                'p1maternal_effect': "", #affects p0
                'p0survival': ""}).lstrip())
        self.model.custom(survival_script)

        self.model.repro(
            population="p1",
            scripts=REPRO_BRYO_DIO_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_BRYO_DIO_P0,
            comment="generates gametes from sporophytes"
            )

        substitution_str = ""
        for i in substitutions:
            substitution_str += "\n  ".join([i.strip(';') + ";\n    "])

        substitution_script = (
            SUBSTITUTION.format(**{'muts': substitution_str,
                'late': LATE_BRYO_DIO}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
            )

    def monoicous(self):
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
            scripts=REPRO_BRYO_MONO_P1,
            comment="generates gametes from sporophytes"
            )

        self.model.repro(
            population="p0",
            scripts=REPRO_BRYO_MONO_P0,
            comment="generates gametes from sporophytes"
            )

        substitution_str = ""
        for i in substitutions:
            substitution_str += "\n  ".join([i.strip(';') + ";\n    "])

        substitution_script = (
            SUBSTITUTION.format(**{'muts': substitution_str,
                'late': LATE_BRYO_MONO}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
            )

    def optimize(
        self, 
        desired_eggs_per_gen = None,
        desired_egg_to_sperm_ratio: float.as_integer_ratio=(1,2)):
        """
        Convenience function to choose optimized params
        """
        gam_female_to_male_ratio = self.gam_female_to_male_ratio
        gam_eggs_per_megaspore = self.gam_eggs_per_megaspore
        gam_sperm_per_microspore = self.gam_sperm_per_microspore
        spo_megaspores_per = self.spo_megaspores_per
        gam_random_death_chance = self.gam_random_death_chance
        fertilization_chance = self.fertilization_chance,
        gam_pop_size = self.gam_pop_size
        
        if not desired_eggs_per_gen:
            desired_eggs_per_gen = (gam_pop_size*gam_female_to_male_ratio
                *spo_megaspores_per*gam_eggs_per_megaspore*
                (1-gam_random_death_chance))

        minimum_sperm_needed = desired_eggs_per_gen/desired_egg_to_sperm_ratio
        minimum_spo_microspores_per = minimum_sperm_needed/gam_sperm_per_microspore

        self.minimum_spo_microspores_per = minimum_spo_microspores_per

        logger.info(f"Optimized `spo_spores_per` is {minimum_spo_microspores_per}.")

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

        mod.reproduction.bryophyte(
            mode='dio',
            spo_pop_size=1000, 
            gam_pop_size=1000,
        )
    print(mod.script)

