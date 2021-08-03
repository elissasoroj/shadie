#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""

from dataclasses import dataclass, field
from shadie.reproduction.base import ReproductionBase
from shadie.reproduction.base_scripts import (
    ACTIVATE, DEACTIVATE, EARLY, SURV,
    SUBSTITUTION, SUB_INNER, REPRO_BRYO_DIO_P1, REPRO_BRYO_DIO_P0,
    REPRO_BRYO_MONO_P1, REPRO_BRYO_MONO_P0, EARLY1_ANGIO,
    REPRO_ANGIO_DIO_P1, REPRO_ANGIO_DIO_P0, REPRO_ANGIO_MONO_P1,
    ANGIO_SURV_P0
)

DTYPES = ("d", "dio", "dioecy", "dioecious", "heterosporous", "dioicous")
MTYPES = ("m", "mono", "monoecy", "monecious", "homosporous", "monoicous")

@dataclass
class AngiospermBase(ReproductionBase):
    lineage: str = field(default="Angiosperm", init=False)
    mode: str
    chromosome: 'shadie.chromosome.ChromosomeBase'

@dataclass
class Angiosperm(AngiospermBase):
    """
    Reproduction mode based on angiosperms; appropriate for gymnosperms
    as well
    """
    diploid_ne: int
    haploid_ne: int
    female_to_male_ratio: float=0.5
    clone_rate: float=0.0
    ovule_count: int=30
    fertilization_rate: float=0.7
    pollen_count: int=100
    pollen_comp: str="F"
    pollen_per_stigma: int=5
    random_death_chance: float=0

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        self.add_initialize_constants()
        self.add_early_haploid_diploid_subpops()        
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
        constants["dK"] = self.diploid_ne
        constants["hK"] = self.haploid_ne
        constants["Death_chance"] = self.random_death_chance
        constants["FtoM"] = self.female_to_male_ratio
        constants["Clone_rate"] = self.clone_rate
        # constants["Clone_num"] = self.clone_number
        # constants["Self_rate"] = self.selfing_rate
        constants["Ovule_count"] = self.ovule_count
        constants["Fertilization_rate"] = self.fertilization_rate
        constants["Pollen_count"] = self.pollen_count
        constants["Pollen_comp"] = self.pollen_comp
        constants["Pollen_count"] = self.clone_rate
        constants["Pollen_per_stigma"] = self.pollen_per_stigma


    def add_early_haploid_diploid_subpops(self):
        """
        add haploid and diploid life stages
        """
        self.model.early(
            time=1,
            scripts= EARLY1_ANGIO,
            comment="define Angiosperm subpops: diploid sporophytes, haploid gametophytes",
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

        early_script = (
            EARLY.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'maternal_effect': "",
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
        """
        fills the script reproduction block with bryophyte-dioicous
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

        early_script = (
            EARLY.format(**{'activate': activate_str, 
                'deactivate': deactivate_str}).lstrip())

        self.model.early(
            time=None, 
            scripts=early_script, 
            comment="alternation of generations",
        )

        survival_script = (
            SURV.format(**{'maternal_effect': "",
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

        mod.reproduction.angiosperm(
            mode='mono',
            chromosome = chrom,
            diploid_ne=1000, 
            haploid_ne=1000,
        )

    print(mod.script)
    #mod.run()
