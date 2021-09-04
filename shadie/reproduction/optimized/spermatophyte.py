#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""
import pyslim
from typing import Union
from dataclasses import dataclass, field
from shadie.reproduction.optimized.base_wf import ReproductionBase
from shadie.reproduction.base_scripts import (
    ACTIVATE, DEACTIVATE, EARLY, SURV,
    SUBSTITUTION, SUB_INNER, REPRO_BRYO_DIO_P1, REPRO_BRYO_DIO_P0,
    REPRO_BRYO_MONO_P1, REPRO_BRYO_MONO_P0, EARLY1_ANGIO,
    REPRO_ANGIO_DIO_P1, REPRO_ANGIO_DIO_P0, REPRO_ANGIO_MONO_P1,
    ANGIO_SURV_P0
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
    gam_pop_size: int
    spo_mutation_rate: Union[None, float] = None
    gam_mutation_rate: float = 0.0
    spo_female_to_male_ratio: float.as_integer_ratio = (1,1)
    spo_clone_rate: float=0.0
    spo_clone_number:int = 1
    flower_ovules_per: int=30
    ovule_fertilization_rate: float=0.7
    pollen_success_rate: float=0.7
    flower_pollen_per: int=100
    pollen_comp: str="F"
    pollen_per_ovule: int=5
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


    def add_initialize_constants(self):
        """
        Add defineConstant calls to init for new variables
        """
        constants = self.model.map["initialize"][0]['constants']
        constants["spo_pop_size"] = self.spo_pop_size
        constants["gam_pop_size"] = self.gam_pop_size
        constants["spo_mutation_rate"] = self.spo_mutation_rate
        constants["gam_mutation_rate"] = self.gam_mutation_rate
        constants["spo_female_to_male_ratio"] = self.spo_female_to_male_ratio
        constants["spo_clone_rate"] = self.spo_clone_rate
        constants["spo_clone_number"] = self.spo_clone_number
        # constants["Self_rate"] = self.selfing_rate
        constants["flower_ovules_per"] = self.flower_ovules_per
        constants["ovule_fertilization_rate"] = self.ovule_fertilization_rate
        constants["pollen_success_rate"] = self.pollen_success_rate
        constants["flower_pollen_per"] = self.flower_pollen_per
        constants["pollen_comp"] = self.pollen_comp
        constants["pollen_per_ovule"] = self.pollen_per_ovule
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
                'p1maternal_effect': "",
                'p0survival': ANGIO_SURV_P0}).lstrip())
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
            SUBSTITUTION.format(**{'inner': substitution_str}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
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
            mode='mono',
            spo_pop_size=1000, 
            gam_pop_size=1000,
        )

    print(mod.script)
    #mod.run()
