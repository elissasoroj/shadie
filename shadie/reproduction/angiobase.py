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
    P0_FITNESS_SCALE_DEFAULT,
    EARLY,
    P1_FITNESS_SCALE_DEFAULT,
    P0_FITNESS_SCALE_DEFAULT,
)
from shadie.reproduction.angio_scripts import (
    REPRO_ANGIO_DIO_P1,
    REPRO_ANGIO_DIO_P0,
    REPRO_ANGIO_MONO_P1,
    REPRO_ANGIO_MONO_P0,
    EARLY1_ANGIO_DIO,
    ANGIO_DIO_FITNESS_SCALE,
    DEFS_ANGIO,
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
    spo_random_death_chance: float
    gam_random_death_chance: float
    spo_maternal_effect: float
    spo_pollen_per: int
    spo_archegonia_per: int
    fertilization_rate: float
    pollen_success_rate: float
    pollen_competition: str
    stigma_pollen_per: int

    def _set_mutation_rates(self):
        """Checks parameters after init."""
        # Set mutation rates for both, or use Model rate for sporo, 0 for gam.
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
                p1_maternal_effect=SPO_MATERNAL_EFFECT_ON_P0)
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

        if self.pollen_competition and self.pollen_competition != "F":
            self.pollen_competition = "T"
        else:
            self.pollen_competition = "F"
            self.stigma_pollen_per = int(1)

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
        self._add_early_script()

    def _define_subpopulations(self):
        """Defines the subpopulations and males/females.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts =[ "p1.individuals.tag=0"])
        else:
            self.model.early(
                time=1,
                scripts=EARLY1_ANGIO_DIO,
                comment="define subpops: p1=diploid sporophytes, p0=haploid gametophytes",
            )

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        early_script = (EARLY.format(
            p0_fitnessScaling = ANGIO_DIO_FITNESS_SCALE,
            p1_fitnessScaling = P1_FITNESS_SCALE_DEFAULT,
            p0activate= self._p0activate_str,
            p0deactivate= self._p0deactivate_str,
            p1activate= self._p1activate_str,
            p1deactivate= self._p1deactivate_str
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="alternation of generations",
        )

    def _add_mode_scripts(self):
        """scripts specific to this organism."""
        #self.model.custom(scripts=FUNCTIONS_ANGIO_MONO, comment = "shadie DEFINITIONS")

        # add reproduction calls
        self.model.repro(
            population="p0",
            scripts=REPRO_ANGIO_DIO_P0,
            idx = "s5",
            comment="generates sporophytes from gametes"
        )
        self.model.repro(
            population="p1",
            scripts=REPRO_ANGIO_DIO_P1,
            idx = "s6",
            comment="generates gametes from sporophytes"
        )
        
        # # add late call
        # substitution_script = (
        #     SUBSTITUTION.format(**{'muts': self._substitution_str,
        #         'late': LATE_ANGIO_DIO}).lstrip())

        # self.model.late(
        #     time=None,
        #     scripts=substitution_script,
        #     comment="fixes mutations in haploid gen"
        #     )


@dataclass
class AngiospermMonoecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""
    mode: str = field(default="monecious", init=False)
    spo_self_rate_per_egg: float
    spo_self_rate: float

    def __post_init__(self):
        #set pollen competition as string
        if self.pollen_competition and self.pollen_competition != "F":
            self.pollen_competition = "T"
        else:
            self.pollen_competition = "F"
            self.stigma_pollen_per = int(1)

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent AngiospermBase class
        self._set_mutation_rates()
        self._add_shared_mode_scripts()

        # methods inherited from parent NonWrightFisher class
        self._add_initialize_constants()
        self._add_alternation_of_generations()
        self._write_trees_file()

        # specific organism functions
        self._define_subpopulations()
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
            p0activate= self._p0activate_str,
            p0deactivate= self._p0deactivate_str,
            p1activate= self._p1activate_str,
            p1deactivate= self._p1deactivate_str
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="alternation of generations",
        )

    def _add_mode_scripts(self):
        """scripts specific to this organism."""
        #self.model.custom(scripts=FUNCTIONS_ANGIO_MONO, comment = "shadie DEFINITIONS")

        self.model.repro(
            population="p0",
            scripts=REPRO_ANGIO_MONO_P0,
            idx = "s5",
            comment="generates sporophytes from gametes"
        )
        self.model.repro(
            population="p1",
            scripts=REPRO_ANGIO_MONO_P1,
            idx = "s6",
            comment="generates gametes from sporophytes"
        )

        # # add late call
        # substitution_script = (
        #     SUBSTITUTION.format(**{'muts': self._substitution_str,
        #         'late': LATE_ANGIO_MONO}).lstrip())

        # self.model.late(
        #     time=None,
        #     scripts=substitution_script,
        #     comment="fixes mutations in haploid gen"
        #     )


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
