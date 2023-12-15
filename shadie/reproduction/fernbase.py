#!/usr/bin/env python

"""
Starting an alternate implementation of Reproduction 
"""
from typing import Union, Optional, Tuple
from dataclasses import dataclass, field
import pyslim
from shadie.reproduction.base import NonWrightFisher
from shadie.reproduction.scripts import (
    EARLY,
    P0_FITNESS_SCALE_DEFAULT,
    P1_FITNESS_SCALE_DEFAULT,
    SPO_CLONES,
    NO_SPO_CLONES,
    SPO_MATERNAL_EFFECT_ON_P0,
    NO_SPO_MATERNAL_EFFECT,
    GAM_CLONES,
    NO_GAM_CLONES,
    GAM_MATERNAL_EFFECT_ON_P1,
    NO_GAM_MATERNAL_EFFECT
)
from shadie.reproduction.fern_scripts import (
    REPRO_PTER_HOMOSPORE_P1, 
    REPRO_PTER_HOMOSPORE_P0,
    REPRO_PTER_HETEROSPORE_P1, 
    REPRO_PTER_HETEROSPORE_P0,
    PTER_FITNESS_SCALE,
    DEFS_PTER_HOMOSPORE,
    DEFS_PTER_HETEROSPORE,
)

from shadie.reproduction.vittaria_scripts import (
    REPRO_PTER_VITTARIA_P0,
    REPRO_PTER_VITTARIA_P1,
    DEFS_PTER_VITTARIA,
    )

DTYPES = ("dioicy", "dioicous", "heterosporous")
MTYPES = ("monoicy", "monoicous", "homosporous")

@dataclass
class PteridophyteBase(NonWrightFisher):
    lineage: str = field(default="Pteridophyte", init=False)
    spo_pop_size: int
    gam_pop_size: int
    spo_mutation_rate: Optional[float]
    gam_mutation_rate: Optional[float]
    spo_clone_rate: float
    spo_clones_per: int
    spo_self_rate: float
    spo_self_rate_per_egg: float
    spo_spores_per: int
    spo_random_death_chance: float
    gam_random_death_chance: float
    spo_maternal_effect: float
    gam_archegonia_per: int
    gam_ceiling: int
    gam_female_to_male_ratio: Tuple[float,float]
    _gens_per_lifecycle: int = field(default=2, init=False)

    def __post_init__(self):
        """Convert tuple ratio to a float."""
        sum_ratio = sum(self.gam_female_to_male_ratio)
        float_ratio = self.gam_female_to_male_ratio[0] / sum_ratio
        self.gam_female_to_male_ratio = float_ratio
        self.model_source = "shadie"
        self.lineage = self.lineage
        self.mode = self.mode
        self.gens_per_lifecycle = self._gens_per_lifecycle

    #TODO?
    #optional (lineage-specific params that correspond to generalized ones)
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

        Adds shadie-defined functions
        """

@dataclass
class PteridophyteHomosporous(PteridophyteBase):
    """Reproduction mode based on homosporoous ferns and lycophytes"""
    mode: str = field(default="homosporous", init=False)
    gam_self_rate: float
    gam_self_rate_per_egg: float
    gam_maternal_effect: float
    gam_clone_rate: float
    gam_clones_per: int

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Pteridophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()
        self._add_first_script()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._set_gametophyte_k()
        self._write_trees_file()
        self._add_initialize_globals()
        self._add_initialize_constants()

        # mode-specific functions
        self._add_mode_scripts()
        self._add_early_script()

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        early_script = (EARLY.format(
            p0_fitnessScaling= P0_FITNESS_SCALE_DEFAULT,
            p1_fitnessScaling= P1_FITNESS_SCALE_DEFAULT,
            gametophyte_clones=GAM_CLONES,
            gam_maternal_effect=GAM_MATERNAL_EFFECT_ON_P1,
            sporophyte_clones=SPO_CLONES,
            spo_maternal_effect=SPO_MATERNAL_EFFECT_ON_P0,
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="alternation of generations",
        )

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to homosporous pteridophyte."""

        self.model.custom(scripts=DEFS_PTER_HOMOSPORE, comment = "shadie DEFINITIONS")

        self.model.repro(
            idx = "s0",
            population="p0",
            scripts=REPRO_PTER_HOMOSPORE_P0,
            comment="generates gametes from sporophytes"
        )
        self.model.repro(
            idx = "s1",
            population="p1",
            scripts=REPRO_PTER_HOMOSPORE_P1,
            comment="generates gametes from sporophytes"
        )

@dataclass
class PteridophyteHeterosporous(PteridophyteBase):
    mode: str = field(default="heterosporous", init=False)

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Bryophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._set_gametophyte_k()
        self._add_initialize_globals()
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
            p0_fitnessScaling=PTER_FITNESS_SCALE,
            p1_fitnessScaling=P1_FITNESS_SCALE_DEFAULT,
            gametophyte_clones=NO_GAM_CLONES,
            gam_maternal_effect=NO_GAM_MATERNAL_EFFECT,
            sporophyte_clones=SPO_CLONES,
            spo_maternal_effect=SPO_MATERNAL_EFFECT_ON_P0,
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="alternation of generations",
        )

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
        self.model.custom(scripts=DEFS_PTER_HETEROSPORE, comment = "shadie DEFINITIONS")

        self.model.repro(
            population="p0",
            idx = "s0",
            scripts=REPRO_PTER_HETEROSPORE_P0,
            comment="generates gametes from sporophytes"
        )
        self.model.repro(
            population="p1",
            idx = "s1",
            scripts=REPRO_PTER_HETEROSPORE_P1,
            comment="generates gametes from sporophytes"
        )
        

@dataclass
class PteridophyteVittaria(PteridophyteBase):
    mode: str = field(default="vittaria", init=False)
    gam_self_rate: float
    gam_self_rate_per_egg: float
    gam_maternal_effect: float
    gam_clone_rate: float
    gam_clones_per: int
    sex: str
    sex_rate: float

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Pteridophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()
        self._add_first_script()

        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._set_gametophyte_k()
        self._add_initialize_globals()
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
            p0_fitnessScaling= P0_FITNESS_SCALE_DEFAULT,
            p1_fitnessScaling= P1_FITNESS_SCALE_DEFAULT,
            gametophyte_clones=GAM_CLONES,
            gam_maternal_effect=GAM_MATERNAL_EFFECT_ON_P1,
            sporophyte_clones=SPO_CLONES,
            spo_maternal_effect=SPO_MATERNAL_EFFECT_ON_P0,
            )
        )

        self.model.early(
            time=None,
            scripts= early_script,
        )

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to Vittaria."""

        self.model.custom(scripts=DEFS_PTER_VITTARIA, comment = "shadie DEFINITIONS")

        self.model.repro(
            population="p0",
            idx = "s0",
            scripts=REPRO_PTER_VITTARIA_P0,
            comment="generates gametes from gametophytes"
        )
        self.model.repro(
            population="p1",
            idx = "s1",
            scripts=REPRO_PTER_VITTARIA_P1,
            comment="generates spores from sporophytes"
        )


if __name__ == "__main__":


    import shadie
    with shadie.Model() as mod:
        
        # define mutation types
        m0 = shadie.mtype(0.5, 'n', (0, 0.4))
        m1 = shadie.mtype(0.5, 'g', (0.8, 0.75), affects_haploid=False)
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

        # mod.reproduction.pteridophyte_vittaria(
        #     spo_pop_size=1000, 
        #     gam_pop_size=1000,
        #     spo_self_rate_per_egg=0.0,
        #     gam_clones_per=10,
        #     gam_clone_rate=0.99,
        #     #spo_spores_per = 100
        # )

        mod.reproduction.pteridophyte_heterosporous(
            spo_pop_size=1000, 
            gam_pop_size=1000,
            spo_self_rate_per_egg=0.0,
            spo_clones_per=2,
            spo_clone_rate=0.02,
            #spo_spores_per = 100
        )


    print(mod.script)
    #print(m1._expr)
    for elem in chrom.elements:
        for mut in elem.mlist:
            print(mut.affects_haploid)


    #mod.run()
