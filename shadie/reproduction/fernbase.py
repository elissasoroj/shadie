#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""
from typing import Union, Optional, Tuple
from dataclasses import dataclass, field
import pyslim
from shadie.reproduction.base import NonWrightFisher
from shadie.reproduction.scripts import (
    SURV,
    GAM_MATERNAL_EFFECT_ON_P1,
    SPO_MATERNAL_EFFECT_ON_P0,
    SUBSTITUTION,
    EARLY,
)
from shadie.reproduction.fern_scripts import (
    REPRO_PTER_HOMOSPORE_P1, 
    REPRO_PTER_HOMOSPORE_P0,
    LATE_PTER_HOMOSPORE,
    REPRO_PTER_HETEROSPORE_P1, 
    REPRO_PTER_HETEROSPORE_P0,
    LATE_PTER_HETEROSPORE,
    PTER_HETERO_FITNESS_SCALE,
    FUNCTIONS_PTER_HOMOSPORE,
    FUNCTIONS_PTER_HETEROSPORE,
    S4_TAG,
)

DTYPES = ("dioicy", "dioicous", "heterosporous")
MTYPES = ("monoicy", "monoicous", "homosporous")

@dataclass
class PteridophyteBase(NonWrightFisher):
    lineage: str = field(default="Bryophyte", init=False)
    spo_pop_size: int
    gam_pop_size: int
    spo_mutation_rate: Optional[float]
    gam_mutation_rate: Optional[float]
    gam_clone_rate: float
    gam_clones_per: int
    spo_clone_rate: float
    spo_clones_per: int
    egg_spo_self_rate: float
    spo_self_chance: float
    spo_random_death_chance: float
    gam_random_death_chance: float
    spo_maternal_effect: float
    gam_archegonia_per: int

    #TODO?
    #optional (lineage-specific params that correspon to generalized ones)
    # cone_megasporangia_per: Optional[int]
    # cone_microsporangia_per: Optional[int]

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
        survival_script = (
                    SURV.format(
                        p0_maternal_effect=SPO_MATERNAL_EFFECT_ON_P0,
                        p1_maternal_effect=GAM_MATERNAL_EFFECT_ON_P1,
                        p0survival="",
                        s4_tag = S4_TAG,
                    ))
        self.model.custom(survival_script, comment="maternal effects and survival")

@dataclass
class PteridophyteHomosporous(PteridophyteBase):
    """Reproduction mode based on homosporoous ferns and lycophytes"""
    mode: str = field(default="homosporous", init=False)
    gam_maternal_effect: float
    gam_self_rate: float
    spo_spores_per: int

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Pteridophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._add_early_script()
        self._add_initialize_constants()
        self._write_trees_file()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
        self.model.custom(scripts=FUNCTIONS_PTER_HOMOSPORE, comment = "shadie DEFINITIONS")

        self.model.repro(
            idx = "s5",
            population="p1",
            scripts=REPRO_PTER_HOMOSPORE_P1,
            comment="generates gametes from sporophytes"
        )
        self.model.repro(
            idx = "s6",
            population="p0",
            scripts=REPRO_PTER_HOMOSPORE_P0,
            comment="generates gametes from sporophytes"
        )
        
        # add late call
        substitution_script = (
            SUBSTITUTION.format(**{'muts': self._substitution_str,
                'late': LATE_PTER_HOMOSPORE}).lstrip())

        self.model.late(
            time=None,
            scripts=substitution_script,
            comment="fixes mutations in haploid gen"
            )

@dataclass
class PteridophyteHeterosporous(PteridophyteBase):
    mode: str = field(default="heterosporous", init=False)
    rs_megasporangia_per: int
    rs_microsporangia_per: int
    megasporangia_megaspores_per: int
    microsporangia_microspores_per: int

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Bryophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._add_initialize_constants()
        self._write_trees_file()

        # mode-specific functions
        self._add_mode_scripts()
        self._add_early_script()

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        early_script = (EARLY.format(
            p0_fitnessScaling= PTER_HETERO_FITNESS_SCALE,
            activate=self._activate_str,
            deactivate=self._deactivate_str
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="alternation of generations",
        )

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
        self.model.custom(scripts=FUNCTIONS_PTER_HETEROSPORE, comment = "shadie DEFINITIONS")

        self.model.repro(
            population="p1",
            idx = "s5",
            scripts=REPRO_PTER_HETEROSPORE_P1,
            comment="generates gametes from sporophytes"
        )
        self.model.repro(
            population="p0",
            idx = "s6",
            scripts=REPRO_PTER_HETEROSPORE_P0,
            comment="generates gametes from sporophytes"
        )
        
        # add late call
        substitution_script = (
            SUBSTITUTION.format(**{'muts': self._substitution_str,
                'late': LATE_PTER_HETEROSPORE}).lstrip())

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
        )


    print(mod.script)
    #mod.run()
