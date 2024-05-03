#!/usr/bin/env python

"""Angiosperm reproduction class is a superclass of NonWrightFisher class.

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
    SPO_CLONES,
    NO_SPO_CLONES,
    GAM_CLONES,
    NO_GAM_CLONES,
    GAM_MATERNAL_EFFECT_ON_P1,
    NO_GAM_MATERNAL_EFFECT,
    SPO_MATERNAL_EFFECT_ON_P0,
    NO_SPO_MATERNAL_EFFECT,
    EARLY,
    P1_FITNESS_SCALE_DEFAULT,
    P0_FITNESS_SCALE_DEFAULT,
)
from shadie.reproduction.angio_scripts import (
    REPRO_ANGIO_DIO_P1,
    REPRO_ANGIO_DIO_P0,
    REPRO_ANGIO_MONO_P1,
    REPRO_ANGIO_MONO_P0,
    FIRST1_ANGIO_MONO,
    FIRST1_ANGIO_DIO,
    ANGIO_DIO_FITNESS_SCALE,
    DEFS_ANGIO_MONO,
    DEFS_ANGIO_DIO,
    POLLEN_COMPETITION,
    NO_POLLEN_COMPETITION,
    NO_EGG_FITNESS,
    EGG_FITNESS_AFFECTS_VIABILITY,
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
    gam_ceiling: int
    _gens_per_lifecycle: int = field(default=2, init=False)

    def __post_init__(self):
        """Add extra params to metadata"""
        self.gens_per_lifecycle = self._gens_per_lifecycle
        self.lineage = self.lineage
        self.mode = self.mode

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

        """


@dataclass
class AngiospermDioecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""
    mode: str = field(default="dioecious", init=False)
    spo_female_to_male_ratio: Tuple[float,float]
    fitness_affects_spo_survival: bool = True
    fitness_affects_spo_reproduction: bool = False
    fitness_affects_gam_survival: bool = True
    fitness_affects_gam_mating: bool = False

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
        self._add_alternation_of_generations()
        self._write_trees_file()
        self._define_subpopulations()

        # specific organism functions
        self._set_gametophyte_k()
        self._add_mode_scripts()
        self._add_early_script()

        #save metadata
        self._add_initialize_globals()
        self._add_initialize_constants()

    def _define_subpopulations(self):
        """Defines the subpopulations and males/females.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts =[ "p1.individuals.tag=0"])
        else:
            self.model.first(
                time=1,
                scripts=FIRST1_ANGIO_DIO,
                comment="define subpops: p1=diploid sporophytes, p0=haploid gametophytes",
            )

    def _add_early_script(self):
        """Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.fitness_affects_gam_survival:
            p0_survival_effects = ANGIO_DIO_FITNESS_SCALE
        else:
            p0_survival_effects = P0_RANDOM_SURVIVAL

        if self.fitness_affects_spo_survival:
            p1_survival_effects = P1_FITNESS_SCALE_DEFAULT
        else:
            p1_survival_effects = P1_RANDOM_SURVIVAL
        early_script = (
            EARLY.format(
                p0_survival_effects = p0_survival_effects,
                p1_survival_effects = p1_survival_effects,
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
        """scripts specific to this organism."""
        self.model.custom(scripts=DEFS_ANGIO_DIO, comment = "shadie DEFINITIONS")

        #add fitness determination of sperm success (or not)
        if self.fitness_affects_gam_mating:
            repro_script_p0 = REPRO_ANGIO_DIO_P0.format(
                pollen_selection=POLLEN_COMPETITION)
        else:
            repro_script_p0 = REPRO_ANGIO_DIO_P0.format(
                pollen_selection=NO_POLLEN_COMPETITION)

        #add fitness determination of spore # (or not)
        if self.fitness_affects_spo_reproduction:
            repro_script_p1 = REPRO_ANGIO_DIO_P1.format(
                egg_selection = EGG_FITNESS_AFFECTS_VIABILITY)

        else: repro_script_p1 = REPRO_ANGIO_DIO_P1.format(
                 egg_selection = NO_EGG_FITNESS)

        # add reproduction calls
        self.model.repro(
            population="p0",
            scripts=repro_script_p0,
            idx="s0",
            comment="generates sporophytes from gametes"
        )
        self.model.repro(
            population="p1",
            scripts=repro_script_p1,
            idx="s1",
            comment="generates gametes from sporophytes"
        )


@dataclass
class AngiospermMonoecious(AngiospermBase):
    """Superclass of Spermatophyte base with dioecious functions."""
    mode: str = field(default="monecious", init=False)
    spo_self_rate_per_egg: float
    spo_self_rate: float
    fitness_affects_spo_survival: bool = True
    fitness_affects_spo_reproduction: bool = False
    fitness_affects_gam_survival: bool = True
    fitness_affects_gam_mating: bool = False

    def __post_init__(self):
        # set pollen competition as string
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
        self._add_alternation_of_generations()
        self._write_trees_file()
        self._set_gametophyte_k()

        # specific organism functions
        self._define_subpopulations()
        self._add_mode_scripts()
        self._add_early_script()

        # add metadata
        self._add_initialize_globals()
        self._add_initialize_constants()

    def _define_subpopulations(self):
        """Defines the subpopulations and males/females.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts =[ "p1.individuals.tag=0"])
        else:
            self.model.first(
                time=1,
                scripts=FIRST1_ANGIO_MONO,
                comment="define subpops: p1=diploid sporophytes, p0=haploid gametophytes",
            )

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.fitness_affects_gam_survival:
            p0_survival_effects = P0_FITNESS_SCALE_DEFAULT
        else:
            p0_survival_effects = P0_RANDOM_SURVIVAL

        if self.fitness_affects_spo_survival:
            p1_survival_effects = P1_FITNESS_SCALE_DEFAULT
        else:
            p1_survival_effects = P1_RANDOM_SURVIVAL
        early_script = (EARLY.format(
            p0_survival_effects = p0_survival_effects,
            p1_survival_effects = p1_survival_effects,
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
        """scripts specific to this organism."""
        self.model.custom(scripts=DEFS_ANGIO_MONO, comment = "shadie DEFINITIONS")
        
        #add fitness determination of sperm success (or not)
        if self.fitness_affects_gam_mating:
            repro_script_p0 = REPRO_ANGIO_MONO_P0.format(
                 pollen_selection=POLLEN_COMPETITION)
        else:
            repro_script_p0 = REPRO_ANGIO_MONO_P0.format(
                 pollen_selection=NO_POLLEN_COMPETITION)

        #add fitness determination of spore # (or not)
        if self.fitness_affects_spo_reproduction:
            repro_script_p1 = REPRO_ANGIO_MONO_P1.format(
                egg_selection = EGG_FITNESS_AFFECTS_VIABILITY)

        else: repro_script_p1 = REPRO_ANGIO_MONO_P1.format(
                egg_selection = NO_EGG_FITNESS)

        self.model.repro(
            population="p0",
            scripts=repro_script_p0,
            idx="s0",
            comment="generates sporophytes from gametes"
        )
        self.model.repro(
            population="p1",
            scripts=repro_script_p1,
            idx="s1",
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
        m0 = shadie.mtype(0.5, 'n', (0, 0.4))
        m1 = shadie.mtype(0.5, 'g', [0.8, 0.75])
        # I suggest we add a checkpoint that calculates the average
        # fitness of mutations input by the user. If fitness is too high
        # the simuulation will lag tremendously.
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

        mod.reproduction.angiosperm_dioecious(
            spo_pop_size=1000,
            gam_pop_size=1000,
            spo_female_to_male_ratio = (1,1)
        )

    print(mod.script)
    # mod.run()
