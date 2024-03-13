#!/usr/bin/env python

"""
Bryophyte reproduction class is a superclass of NonWrightFisher class.

Class inheritance structure
---------------------------
ReproductionBase
    WrightFisher
    NonWrightFisher
        BrypophyteBase
            BryophyteDioicous
            BryophyteMonoicous
"""

from typing import Tuple, Optional
from dataclasses import dataclass, field
from shadie.reproduction.base import NonWrightFisher
from shadie.reproduction.scripts import (
    SPO_CLONES,
    NO_SPO_CLONES,
    GAM_CLONES,
    NO_GAM_CLONES,
    GAM_MATERNAL_EFFECT_ON_P1,
    NO_GAM_MATERNAL_EFFECT,
    SPO_MATERNAL_EFFECT_ON_P0,
    NO_SPO_MATERNAL_EFFECT,
    P0_FITNESS_SCALE_DEFAULT,
    EARLY,
    P0_FITNESS_SCALE_DEFAULT,
    P1_FITNESS_SCALE_DEFAULT,
)
from shadie.reproduction.bryo_scripts import (
    REPRO_BRYO_DIO_P1,
    REPRO_BRYO_DIO_P0,
    REPRO_BRYO_MONO_P1,
    REPRO_BRYO_MONO_P0,
    DEFS_BRYO_MONO,
    DEFS_BRYO_DIO,
)

DTYPES = ("dioicy", "dioicous", "heterosporous")
MTYPES = ("monoicy", "monoicous", "homosporous")

@dataclass
class BryophyteBase(NonWrightFisher):
    """Reproduction mode based on mosses, hornworts, and liverworts."""
    lineage: str = field(default="Bryophyte", init=False)
    spo_pop_size: int
    gam_pop_size: int
    spo_mutation_rate: Optional[float]
    gam_mutation_rate: Optional[float]
    gam_clone_rate: float
    gam_clones_per: int
    #gam_self_rate: float
    spo_self_rate_per_egg: float
    #spo_self_rate: float
    spo_random_death_chance: float
    gam_random_death_chance: float
    spo_spores_per: int
    gam_maternal_effect: float
    gam_archegonia_per: int
    gam_ceiling: int
    _gens_per_lifecycle: int = field(default=2, init=False)

    def __post_init__(self):
        """Add extra params to metadata"""
        self.model_source = "shadie"
        self.lineage = self.lineage
        self.mode = self.mode
        self.gens_per_lifecycle = self._gens_per_lifecycle

        """Convert tuple ratio to a float."""
        sum_ratio = sum(self.gam_female_to_male_ratio)
        float_ratio = self.gam_female_to_male_ratio[0] / sum_ratio
        self.gam_female_to_male_ratio = float_ratio

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
        """

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This will be overridden by any callbacks of the same name in subclasses
        """
        early_script = (EARLY.format(
            p0_fitnessScaling= P0_FITNESS_SCALE_DEFAULT,
            p1_fitnessScaling= P1_FITNESS_SCALE_DEFAULT,
            gametophyte_clones=GAM_CLONES,
            gam_maternal_effect=GAM_MATERNAL_EFFECT_ON_P1,
            sporophyte_clones=NO_SPO_CLONES,
            spo_maternal_effect=NO_SPO_MATERNAL_EFFECT,
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="events after reproduction",
        )


@dataclass
class BryophyteDioicous(BryophyteBase):
    mode: str = field(default="dioicous", init=False)
    gam_female_to_male_ratio: Tuple[float,float]

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Bryophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()
        self._add_early_script()

        # methods inherited from parent NonWrightFisher class
        self._add_first_script()
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._set_gametophyte_k()
        self._add_initialize_globals()
        self._add_initialize_constants()
        self._write_trees_file()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
        self.model.custom(scripts=DEFS_BRYO_DIO, comment = "shadie DEFINITIONS")
        self.model.repro(
            population="p0",
            scripts=REPRO_BRYO_DIO_P0,
            idx = "s0",
            comment="generates sporophytes from gametes"
        )
        self.model.repro(
            population="p1",
            scripts=REPRO_BRYO_DIO_P1,
            idx = "s1",
            comment="generates gametes from sporophytes"
        )


@dataclass
class BryophyteMonoicous(BryophyteBase):
    mode: str = field(default="monoicous", init=False)
    gam_self_rate_per_egg: float
    gam_female_to_male_ratio: Tuple[float,float]

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent Bryophyte class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()
        self._add_early_script()

        # methods inherited from parent NonWrightFisher class
        self._add_first_script()
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._set_gametophyte_k()
        self._add_initialize_globals()
        self._add_initialize_constants()
        self._write_trees_file()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """fills the model.map block with bryophyte-monoicous scripts."""
        # add reproduction scripts
        self.model.custom(scripts=DEFS_BRYO_MONO, comment = "shadie DEFINITIONS")
        self.model.repro(
            population="p0",
            scripts=REPRO_BRYO_MONO_P0,
            idx = "s0",
            comment="generates sporophytes from gametes"
        )
        self.model.repro(
            population="p1",
            scripts=REPRO_BRYO_MONO_P1,
            idx="s1",
            comment="generates gametes from sporophytes"
        )

if __name__ == "__main__":

    import shadie

    # define mutation types
    m0 = shadie.mtype(0.5, 'n', [0, 0.4])
    m1 = shadie.mtype(0.5, 'g', [0.8, 0.75])
    m2 = shadie.mtype(0.5, 'g', [0.8, 0.75], affects_diploid = False)
    m3 = shadie.mtype(0.5, 'n', [0, 0.4], affects_haploid = False)

    # define elements types
    e0 = shadie.etype([m0, m2], [1, 2])
    e1 = shadie.etype([m3], [1])

    # design chromosome of elements
    chrom = shadie.chromosome.random(
        genome_size=20000,
        noncds=e0,
        intron=e0,
        exon=e1,
    )

    with shadie.Model() as mod:
        mod.initialize(chromosome=chrom, sim_time=50, file_out="/tmp/test.trees")
        mod.reproduction.bryophyte_monoicous(
            spo_pop_size=100,
            gam_pop_size=100,
            gam_self_rate_per_egg=0.8,
        )
    print(mod.script)
    #mod.write("/tmp/slim.slim")
    #mod.run(seed=123)
