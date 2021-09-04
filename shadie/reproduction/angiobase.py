#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""

from typing import Tuple, Optional
from dataclasses import dataclass, field
from shadie.reproduction.base import NonWrightFisher
from shadie.reproduction.base_scripts import SURV
from shadie.reproduction.angio_scripts import (
    EARLY1_ANGIO,
    REPRO_ANGIO_DIO_P1, 
    REPRO_ANGIO_DIO_P0, 
    REPRO_ANGIO_MONO_P1,
    ANGIO_SURV_P0
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
    spo_flowers_per_individual: int
    spo_ovules_per_flower: int
    spo_ovule_success_rate: float
    spo_pollen_per_flower: int
    spo_pollen_success_rate: float
    spo_pistils_per_flower: int
    spo_pistils_pollen_comp: bool
    spo_clone_rate: float
    spo_clone_number: int
    spo_self_rate: float
    spo_maternal_effect: float
    spo_random_death_chance: float
    gam_random_death_chance: float

    def __post_init__(self):
        """Checks parameters after init."""
        # Set mutation rates for both, or use Model rate / 2 for both.
        if self.spo_mutation_rate or self.gam_mutation_rate:
            require_spo = self.spo_mutation_rate is not None
            require_gam = self.gam_mutation_rate is not None
            assert require_gam and require_spo, (
                "You must define a mutation rate for both sporophyte "
                "and gametophyte generations.")
        else:
            self.spo_mutation_rate = 0.5 * self.model._mutation_rate
            self.gam_mutation_rate = 0.5 * self.model._mutation_rate

    def _add_shared_mode_scripts(self):
        """Survival and maternal effects shared by spermatophytes."""
        survival_script = SURV.format(
            p0maternal_effect="",
            p1maternal_effect="",
            p0survival=ANGIO_SURV_P0,
        )
        self.model.custom(survival_script) 


@dataclass
class AngiospermDioecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""
    spo_female_to_male_ratio: Tuple[float,float]

    def __post_init__(self):
        """set ratio as a float."""
        ratio_sum = sum(self.spo_female_to_male_ratio)
        self.spo_female_to_male_ratio = (
            self.spo_female_to_male_ratio[0] / ratio_sum)

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent NonWrightFisher class
        self._add_initialize_constants()
        self._add_alternation_of_generations()
        self._write_trees_file()

        # methods inherited from AngiospermBase
        self._add_shared_mode_scripts()

        # specific organism functions
        self._define_subpopulations()
        self._add_mode_scripts()

    def _define_subpopulations(self):
        """Defines the subpopulations and males/females.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.model._file_in:
            self.model.read_from_file()
        else:
            self.model.early(
                time=1,
                scripts=EARLY1_ANGIO,
                comment="define subpops and initial sex ratio",
            )

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


@dataclass
class AngiospermMonoecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_initialize_constants()
        self._add_alternation_of_generations()
        self._write_trees_file()

        # methods inherited from AngiospermBase
        self._add_shared_mode_scripts()

        # specific organism functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """scripts specific to this organism."""
        self.model.repro(
            population="p1",
            scripts=REPRO_ANGIO_MONO_P1,
            comment="generates gametes from sporophytes"
        )
        # TODO: seems like this should say MONO not DIO, right? 
        # or are they the same in this case?
        self.model.repro(
            population="p0",
            scripts=REPRO_ANGIO_DIO_P0,
            comment="generates gametes from sporophytes"
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
            spo_popsize=1000, 
            gam_popsize=1000,
        )

    print(mod.script)
    #mod.run()
