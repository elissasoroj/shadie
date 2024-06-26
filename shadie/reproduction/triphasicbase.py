#!/usr/bin/env python

"""Classes for Polysiphonia (red algae, triphasic) reproduction.

Parameters
----------

"""
from typing import Optional, Tuple
from dataclasses import dataclass, field
from shadie.reproduction.base import NonWrightFisher
from shadie.reproduction.triphasic_scripts import (
    FIRST,
    EARLY,
    P1_FITNESS_SCALE_DEFAULT,
    P2_FITNESS_SCALE_DEFAULT,
    TSPO_CLONES,
    NO_TSPO_CLONES,
    GAM_CLONES,
    NO_GAM_CLONES,
    GAM_MATERNAL_EFFECT_ON_P0,
    NO_GAM_MATERNAL_EFFECT,
    FITNESS_AFFECTS_TSPO_REPRODUCTION,
    CONSTANT_TSPORES,
    FITNESS_AFFECTS_SPERM_SUCCESS,
    RANDOM_MATING,
)
from shadie.reproduction.polysiphonia_scripts import (
    DEFS_PSIPH,
    REPRO_PSIPH_P2,
    REPRO_PSIPH_P1,
    REPRO_PSIPH_P0,
    CONSTANT_CARPOGONIA,
    MATERNAL_FITNESS_AFFECTS_CARPOGONIA_NUM,
    CARPOGONIUM_FITNESS_AFFECTS_CARPOSPORE_NUM,
    CONSTANT_CSPORES,
)


@dataclass
class TriphasicBase(NonWrightFisher):
    lineage: str = field(default="Triphasic", init=False)
    gam_pop_size: int
    tspo_pop_size: int
    gam_mutation_rate: Optional[float]
    cspo_mutation_rate: Optional[float]
    tspo_mutation_rate: Optional[float]
    gam_clone_rate: float
    gam_clones_per: int
    tspo_clone_rate: float
    tspo_clones_per: int
    tspo_self_rate_per_egg: float
    tspo_max_spores_per: int
    cspo_max_spores_per: int
    gam_max_carpogonia_per: int
    gam_female_to_male_ratio: Tuple[float, float]
    gam_random_death_chance: float
    cspo_random_death_chance: float
    tspo_random_death_chance: float
    gam_maternal_effect: float
    gam_ceiling: int
    fitness_affects_tspo_survival: bool = True
    fitness_affects_gam_survival: bool = True
    fitness_affects_tspo_reproduction: bool = False
    fitness_affects_sperm_success: bool = False
    fitness_affects_egg_num: bool = False
    cspo_recombination: str = "F"
    _gens_per_lifecycle: int = field(default=3, init=False)

    def __post_init__(self):
        """Convert tuple ratio to a float."""
        sum_ratio = sum(self.gam_female_to_male_ratio)
        float_ratio = self.gam_female_to_male_ratio[0] / sum_ratio

        # save params for metadata output
        self.gam_female_to_male_ratio = float_ratio
        self.model_source = "shadie"
        self.lineage = self.lineage
        self.mode = self.mode
        self.gens_per_lifecycle = self._gens_per_lifecycle

    # TODO?
    # optional (lineage-specific params that correspond to generalized ones)
    # cone_megasporangia_per: Optional[int]
    # cone_microsporangia_per: Optional[int]

    def _set_mutation_rates(self):
        """Checks parameters after init."""
        # Set mutation rates for both, or use Model rate / 2 for both.

        if self.tspo_mutation_rate or self.gam_mutation_rate or self.cspo_mutation_rate:
            require_tspo = self.tspo_mutation_rate is not None
            require_gam = self.gam_mutation_rate is not None
            require_cspo = self.cspo_mutation_rate is not None
            assert require_gam and require_tspo and require_cspo, (
                "You must define a mutation rate for both sporophyte "
                "and gametophyte generations.")
        else:
            self.tspo_mutation_rate = self.model.metadata['mutation_rate']/3
            self.cspo_mutation_rate = self.model.metadata['mutation_rate']/3
            self.gam_mutation_rate = self.model.metadata['mutation_rate']/3

    def _define_subpopulations(self):
        """add haploid and diploid life stages as subpopulations.
        """
        # load a trees file that already has p1 and p2 pops
        if self.model.metadata['file_in']:
            # set p2 to tag=3 (sporophytes) and p1 to tag=[0,1] (gametophytes)
            self.model._read_from_file(tag_scripts=[
                "p2.individuals.tag = 3;",
                # fix this so that maternal fitness is read in from a vector?
                "p2.individuals.setValue('maternal_fitness', 1.0);",
                "tags = (runif(p1.individualCount) < GAM_FEMALE_TO_MALE_RATIO); ",
                "p1.individuals.tagL0 = tags;",
            ])

        # create new p1 and p2 populations
        else:
            self.model.first(
                time=1,
                scripts=[
                    "sim.addSubpop('p2', TSPO_POP_SIZE)",
                    "sim.addSubpop('p1', 0)",
                    "sim.addSubpop('p0', 0)",
                    "p2.individuals.tag = 3",
                    "p2.individuals.setValue('maternal_fitness', 1.0);",
                ],
                comment="define subpops: p2=diploid sporophytes, p1=haploid gametophytes",
            )

    def _add_first_script(self):
        """Defines the first() callbacks for each gen.

        This is overridden by callbacks of the same name in subclasses
        """
        self.model.first(
            time=None,
            scripts=FIRST,
            comment="alternation of generations",
        )

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        if self.fitness_affects_gam_survival:
            p1_survival_effects = P1_FITNESS_SCALE_DEFAULT
        else:
            p1_survival_effects = P1_RANDOM_SURVIVAL
            self.model.survival(
                comment = "fitness doesn't affect gametophyte survival;",
                population = "p1",
                scripts = "return T;"
            )

        if self.fitness_affects_tspo_survival:
            p2_survival_effects = P2_FITNESS_SCALE_DEFAULT
        else:
            p2_survival_effects = P2_RANDOM_SURVIVAL
            self.model.survival(
                comment = "fitness doesn't affect tetrasporophyte survival;",
                population = "p2",
                scripts = "return T;"
            )

        early_script = (
            EARLY.format(
                p1_survival_effects= p1_survival_effects,
                p2_survival_effects= p2_survival_effects,
                gametophyte_clones=GAM_CLONES,
                gam_maternal_effect=GAM_MATERNAL_EFFECT_ON_P0,
                tetrasporophyte_clones=TSPO_CLONES,
            )
        )
        self.model.early(
            time=None,
            scripts=early_script,
        )

@dataclass
class TriphasicPolysiphonia(TriphasicBase):
    """Reproduction mode based on homosporoous ferns and lycophytes"""
    mode: str = field(default="Polysiphonia", init=False)

    def run(self):
        """Fill self.model.map with SLiM script snippets."""
        # methods inherited from parent TriphasicBase class
        self._set_mutation_rates()
        self._add_early_script()
        self._add_first_script()
        self._define_subpopulations()

        # methods inherited from parent NonWrightFisher class
        self._add_alternation_of_generations()
        self._set_gametophyte_k()
        self._write_trees_file()
        self._add_initialize_globals()
        self._add_initialize_constants()

        # mode-specific functions
        self._add_mode_scripts()

    def _add_mode_scripts(self):
        """Add reproduction scripts unique to homosporous pteridophyte."""
        self.model.custom(scripts=DEFS_PSIPH, comment="shadie DEFINITIONS")
        
        #add fitness determination of sperm success (or not)
        if self.fitness_affects_sperm_success:
            if self.fitness_affects_egg_num:
                repro_script_p1 = REPRO_PSIPH_P1.format(
                    carpogonia_determination=MATERNAL_FITNESS_AFFECTS_CARPOGONIA_NUM,
                    sperm_sampling=FITNESS_AFFECTS_SPERM_SUCCESS)
            else:
                repro_script_p1 = REPRO_PSIPH_P1.format(
                    carpogonia_determination=CONSTANT_CARPOGONIA,
                    sperm_sampling=FITNESS_AFFECTS_SPERM_SUCCESS)

        else:
            if self.fitness_affects_egg_num:
                repro_script_p1 = REPRO_PSIPH_P1.format(
                    carpogonia_determination=MATERNAL_FITNESS_AFFECTS_CARPOGONIA_NUM,
                    sperm_sampling=RANDOM_MATING)
            else:
                repro_script_p1 = REPRO_PSIPH_P1.format(
                    carpogonia_determination=CONSTANT_CARPOGONIA,
                    sperm_sampling=RANDOM_MATING)

        #add fitness determination of spore # (or not)
        if self.fitness_affects_tspo_reproduction:
            repro_script_p2 = REPRO_PSIPH_P2.format(
                tetraspore_determination=FITNESS_AFFECTS_TSPO_REPRODUCTION)

        else: repro_script_p2 = REPRO_PSIPH_P2.format(
                tetraspore_determination=CONSTANT_TSPORES)

        #add fitness determination of spore # (or not)
        if self.fitness_affects_tspo_reproduction:
            repro_script_p0 = REPRO_PSIPH_P0.format(
                carpospore_determination=CARPOGONIUM_FITNESS_AFFECTS_CARPOSPORE_NUM)

        else: repro_script_p0 = REPRO_PSIPH_P0.format(
                carpospore_determination=CONSTANT_CSPORES)

        self.model.repro(
            population="p1",
            scripts=repro_script_p1,
            idx = "s1",
            comment="generates carposporophytes from gametes"
        )

        self.model.repro(
            population="p2",
            scripts=repro_script_p2,
            idx = "s2",
            comment="generates tetraspores from tetrasporophytes"
        )

        self.model.repro(
            population="p0",
            scripts=repro_script_p0,
            idx = "s0",
            comment="generates carpospores from carposporophytes"
        )

        self.model.survival(
            comment = "fitness doesn't affect carposporophyte survival;",
            population = "p0",
            scripts = "return T;"
        )

if __name__ == "__main__":


    import shadie
    with shadie.Model() as mod:

        # define mutation types
        m0 = shadie.mtype(0.5, 'n', (0, 0.4))
        m1 = shadie.mtype(0.5, 'g', (0.8, 0.75), affects_haploid=False)
        # I suggest we add a checkpoint that calculates the average
        # fitness of mutations input by the user. If fitness is too high
        # the simuulation will lag tremendously.
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
        mod.initialize(
            chromosome=chrom,
            sim_time=20,
        )

        mod.reproduction.triphasic_polysiphonia(
            tspo_pop_size=500,
            gam_pop_size=1_000,
            tspo_max_spores_per=100,
            gam_max_carpogonia_per=10,
            gam_female_to_male_ratio=(1, 1),
            tspo_clone_rate=0,
            tspo_clones_per=0,
            gam_clone_rate=0,
            gam_clones_per=0,
            tspo_self_rate_per_egg=0,
            tspo_random_death_chance=0,
            gam_random_death_chance=0,
            gam_maternal_effect=0.0,
            gam_ceiling=3_000,
        )

        # mod.reproduction.pteridophyte_heterosporous(
        #     spo_pop_size=1000,
        #     gam_pop_size=1000,
        #     spo_self_rate_per_egg=0.0,
        #     spo_clones_per=2,
        #     spo_clone_rate=0.02,
        #     #spo_spores_per = 100
        # )

    print(mod.script)
    # print(m1._expr)
    for elem in chrom.elements:
        for mut in elem.mlist:
            print(mut.affects_haploid)
    # mod.run()
