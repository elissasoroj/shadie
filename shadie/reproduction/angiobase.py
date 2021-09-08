#!/usr/bin/env python

"""
Angiosperm reproduction class is a superclass of NonWrightFisher class.

Class inheritance structure
---------------------------
ReproductionBase
    WrightFisher
    NonWrightFisher
        AngiospermBase
            AngiospermDioecious
            AngiospermMonecious
"""

from typing import Tuple, Optional, Union
from dataclasses import dataclass, field
from shadie.reproduction.base import NonWrightFisher
from shadie.reproduction.scripts import (
    SURV,
    SPO_MATERNAL_EFFECT_ON_P0,
    SUBSTITUTION,
)
from shadie.reproduction.angio_scripts_opt import (
    REPRO_ANGIO_DIO_P1,
    REPRO_ANGIO_DIO_P0,
    LATE_ANGIO_DIO,
    REPRO_ANGIO_MONO_P1,
    REPRO_ANGIO_MONO_P0,
    LATE_ANGIO_MONO,
    EARLY1_ANGIO,
    ANGIO_P0_SURV,
)

DTYPES = ("d", "dio", "dioecy", "dioecious",)
MTYPES = ("m", "mono", "monoecy", "monecious",)

@dataclass
class AngiospermBase(NonWrightFisher):
    lineage: str = field(default="Angiosperm", init=False)
    spo_pop_size: int
    gam_pop_size: int
    spo_mutation_rate: Optional[float]
    gam_mutation_rate: Optional[float]
    spo_clone_rate: float
    spo_clones_per: int
    spo_self_rate: float
    spo_random_death_chance: float
    gam_random_death_chance: float
    spo_maternal_effect: float
    spo_flowers_per: int
    flower_ovules_per: int
    flower_anthers_per: int
    anther_pollen_per: int
    spo_ovule_success_rate: float
    spo_pollen_success_rate: float
    pollen_comp: Union[bool, str]
    pollen_comp_stigma_pollen_per: int

    def _set_mutation_rates(self):
        """Checks parameters after init."""
        # Set mutation rates for both, or use Model rate / 2 for both.
        if self.spo_mutation_rate or self.gam_mutation_rate:
            require_spo = self.spo_mutation_rate is not None
            require_gam = self.gam_mutation_rate is not None
            assert require_gam and require_spo, (
                "You must define a mutation rate for both sporophyte "
                "and gametophyte generations.")
        else:
            self.spo_mutation_rate = 0.5 * self.model.metadata['mutation_rate']
            self.gam_mutation_rate = 0.5 * self.model.metadata['mutation_rate']

    def _add_shared_mode_scripts(self):
        """Adds scripts shared by homosp and heterosp superclasses.

        Adds a survival script to define the random_chance_of_death,
        maternal effects, and survival=0 for alternation of generations.
        """
        survival_script = (
            SURV.format(
                p0_maternal_effect="",
                p1_maternal_effect=SPO_MATERNAL_EFFECT_ON_P0,
                p0survival=ANGIO_P0_SURV,
            )
        )
        self.model.custom(survival_script, comment="maternal effects and survival")


@dataclass
class AngiospermDioecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""
    mode: str = field(default="dioecious", init=False)
    spo_female_to_male_ratio: Tuple[float,float]

    def __post_init__(self):
        """set ratio as a float."""
        ratio_sum = sum(self.spo_female_to_male_ratio)
        self.spo_female_to_male_ratio = (
            self.spo_female_to_male_ratio[0] / ratio_sum)

        if self.pollen_comp and self.pollen_comp != "F":
            self.pollen_comp = "T"
        else:
            self.pollen_comp = "F"
            self.pollen_comp_stigma_pollen_per = int(1)

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from AngiospermBase
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._add_initialize_constants()
        self._add_alternation_of_generations()
        self._write_trees_file()

        # specific organism functions
        self._define_subpopulations()
        self._add_mode_scripts()

    def _define_subpopulations(self):
        """Defines the subpopulations and males/females.
        This overrides the NonWrightFisher class function of same name.
        """

    def _add_mode_scripts(self):
        """scripts specific to this organism."""
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
        
        # add late call
        substitution_script = (
            SUBSTITUTION.format(**{'muts': self._substitution_str,
                'late': LATE_ANGIO_DIO}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
            )


@dataclass
class AngiospermMonoecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""
    mode: str = field(default="monecious", init=False)

    def __post_init__(self):
        #set pollen competition as string
        if self.pollen_comp and self.pollen_comp != "F":
            self.pollen_comp = "T"
        else:
            self.pollen_comp = "F"
            self.pollen_comp_stigma_pollen_per = int(1)

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent NonWrightFisher class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._add_initialize_constants()
        self._add_alternation_of_generations()
        self._write_trees_file()

        # specific organism functions
        self._define_subpopulations()
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """scripts specific to this organism."""
        self.model.repro(
            population="p1",
            scripts=REPRO_ANGIO_MONO_P1,
            comment="generates gametes from sporophytes"
        )
        self.model.repro(
            population="p0",
            scripts=REPRO_ANGIO_MONO_P0,
            comment="generates gametes from sporophytes"
        )
        
        # add late call
        substitution_script = (
            SUBSTITUTION.format(**{'muts': self._substitution_str,
                'late': LATE_ANGIO_MONO}).lstrip())

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

        mod.reproduction.angiosperm_monoecious(
            spo_pop_size=1000, 
            gam_pop_size=1000,
        )

    print(mod.script)
    #mod.run()
