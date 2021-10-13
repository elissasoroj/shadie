#!/usr/bin/env python

"""Classes for Fern reproduction.

The parent PteridophyteBase class has shared parameters and functions
for classes PteridophyteHomosporous and PteridophyteHeterosporous. The
parameters for these classes are defined in the factory functions 
in :mod:`reproduction.optimized.pteridophytes`.

Parameters
----------

"""
from typing import Optional
from dataclasses import dataclass, field
from shadie.reproduction.base import NonWrightFisher
# from shadie.reproduction.scripts import (
    # GAM_MATERNAL_EFFECT_ON_P1,
    # SPO_MATERNAL_EFFECT_ON_P0,
    # SUBSTITUTION,
    # EARLY,
# )
from shadie.reproduction.fern_scripts import (
    REPRO_PTER_HOMOSPORE_P1, 
    REPRO_PTER_HOMOSPORE_P0,
    LATE_PTER_HOMOSPORE,
    FUNCTIONS_PTER_HOMOSPORE,   
    # PTER_HETERO_FITNESS_SCALE,
)

from shadie.reproduction.fern_scripts2 import (
    LATE_PTER_HETERO,
    EARLY_PTER_HETERO,
    FUNCTIONS_PTER_HETERO,
    SURVIVAL_PTER_HETERO,
    REPRO_PTER_HETEROSPORE_P1, 
    REPRO_PTER_HETEROSPORE_P0,
)


@dataclass
class PteridophyteBase(NonWrightFisher):
    lineage: str = field(default="Bryophyte", init=False)
    spo_pop_size: int
    gam_pop_size: int
    spo_mutation_rate: Optional[float]
    gam_mutation_rate: Optional[float]
    gam_clone_rate: float   # TODO: maybe replace rate/num with a single poisson draw
    gam_clones_per: int
    spo_clone_rate: float
    spo_clones_per: int
    # spo_self_rate: float
    spo_self_rate_per_egg: float  # egg_spo_self_rate: float
    spo_self_rate_per_ind: float  # spo_self_chance: float
    spo_random_death_chance: float
    gam_random_death_chance: float
    spo_maternal_effect: float
    gam_archegonia_per: int
    # gam_antheridia_per: int  # not modeled b/c it is very large.

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

        Adds shadie-defined functions and a survival script to define 
        the random_chance_of_death, maternal effects, and survival=0 for 
        alternation of generations.
        """
        self.model.custom(
            scripts=SURVIVAL_PTER_HETERO,
            comment="alternate to other generation, random death and maternal effects.",
        )


@dataclass
class PteridophyteHomosporous(PteridophyteBase):
    """Reproduction mode based on homosporoous ferns and lycophytes"""
    mode: str = field(default="homosporous", init=False)
    spo_spores_per: int
    gam_maternal_effect: float
    gam_self_rate: float

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Pteridophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_early_and_late(...)
        self._add_initialize_constants()
        self._write_trees_file()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
        self.model.custom(
            scripts=FUNCTIONS_PTER_HOMOSPORE, 
            comment="shadie DEFINITIONS",
        )
        self.model.repro(
            idx="s5",
            population="p1",
            scripts=REPRO_PTER_HOMOSPORE_P1,
            comment="generates gametes from sporophytes"
        )
        self.model.repro(
            idx="s6",
            population="p0",
            scripts=REPRO_PTER_HOMOSPORE_P0,
            comment="generates gametes from sporophytes"
        )


@dataclass
class PteridophyteHeterosporous(PteridophyteBase):
    mode: str = field(default="heterosporous", init=False)
    # rs_megasporangia_per: int
    # rs_microsporangia_per: int
    # megasporangia_megaspores_per: int
    # microsporangia_microspores_per: int
    #spo_female_to_male_ratio: float
    spo_microspores_per: int
    spo_megaspores_per: int

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Bryophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_early_and_late(EARLY_PTER_HETERO, LATE_PTER_HETERO)
        self._add_initialize_constants()
        self._write_trees_file()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
        self.model.custom(
            scripts=FUNCTIONS_PTER_HETERO, 
            comment="shadie DEFINITIONS",
        )

        self.model.repro(
            population="p1",
            idx="s5",
            scripts=REPRO_PTER_HETEROSPORE_P1,
            comment="reproduction in p1: eggs or sperm are created."
        )

        self.model.repro(
            population="p0",
            idx="s6",
            scripts=REPRO_PTER_HETEROSPORE_P0,
            comment="reproduction in p0: eggs are fertilized by sperm."
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
        e0 = shadie.EXON
        e1 = shadie.INTRON
        
        # design chromosome of elements
        chrom = shadie.chromosome.random(
            genome_size=20000,
            noncds=e0,
            intron=e0,
            exon=e1,
        )

        # init the model
        mod.initialize(chromosome=chrom)

        mod.reproduction.pteridophyte_homosporous(
            spo_pop_size=1000, 
            gam_pop_size=1000,
            spo_spores_per=100
        )


    print(mod.script)
    #mod.run()
