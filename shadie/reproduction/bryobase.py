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
from shadie.reproduction.base_scripts import (
    SURV,
    MATERNAL_EFFECT,
)
from shadie.reproduction.bryo_scripts import (
    REPRO_BRYO_DIO_P1,
    REPRO_BRYO_DIO_P0,
    REPRO_BRYO_MONO_P1,
    REPRO_BRYO_MONO_P0,
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
    spo_spores_per: int
    gam_sporophytes_per: int
    gam_clone_rate: float
    gam_clone_number: int
    gam_maternal_effect: float
    spo_self_rate: float
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
        """Fills model.map scripts applying to heterosp or homosp bryos.

        This will define survival callback functions s5-s8 which are
        turned on or off depending on model parameters.
        """
        # add a survival script to define the random_chance_of_death,
        # maternal effects, and survival=0 for alternation of generations.
        survival_script = (
            SURV.format(
                p0maternal_effect="",
                p1maternal_effect=MATERNAL_EFFECT,
                p0survival="return NULL;"
            )
        )
        self.model.custom(survival_script, comment="maternal effects and survival")


@dataclass
class BryophyteDioicous(BryophyteBase):
    mode: str = field(default="heterosporous", init=False)
    gam_female_to_male_ratio: Tuple[float,float]=(2, 1)
    spo_megaspores_per: int=1
    spo_microspores_per: int=1

    def __post_init__(self):
        """Convert tuple ratio to a float."""
        sum_ratio = sum(self.gam_female_to_male_ratio)
        float_ratio = self.gam_female_to_male_ratio[0] / sum_ratio
        self.gam_female_to_male_ratio = float_ratio

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()
        self._add_initialize_constants()
        self._write_trees_file()

        # methods inherited from parent Bryophyte class
        self._add_shared_mode_scripts()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to heterosporous bryo."""
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


@dataclass
class BryophyteMonoicous(BryophyteBase):
    mode: str = field(default="homosporous", init=False)
    gam_self_rate: float=0.2

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent NonWrightFisher class
        self._define_subpopulations()
        self._add_alternation_of_generations()        
        self._add_initialize_constants()
        self._write_trees_file()

        # methods inherited from parent Bryophyte class
        self._add_shared_mode_scripts()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """fills the model.map block with bryophyte-monoicous scripts."""
        # add reproduction scripts
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



if __name__ == "__main__":

    import shadie

    # define mutation types
    m0 = shadie.mtype(0.5, 'n', 0, 0.4)
    m1 = shadie.mtype(0.5, 'g', 0.8, 0.75)

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

    with shadie.Model() as mod:
        mod.initialize(chromosome=chrom, sim_time=50, file_out="/tmp/test.trees")
        mod.reproduction.bryophyte_dioicous(
            spo_pop_size=100,
            gam_pop_size=100,
        )
    mod.write("/tmp/slim.slim")
    # mod.run(binary="/usr/local/bin/slim", seed=123)
