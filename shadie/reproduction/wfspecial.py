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
from shadie.reproduction.base import ReproductionBase
from shadie.reproduction.base import WrightFisher
from shadie.reproduction.scripts import (
    GAM_MATERNAL_EFFECT_ON_P1,
    P0_FITNESS_SCALE_DEFAULT,
    #EARLY_WITH_GAM_K,
    EARLY,
    HAP_MUT_FITNESS,
    DIP_MUT_FITNESS
)
from shadie.reproduction.specialWF_scripts import (
    REPRO_HAPLOID_WF,
    REPRO_CLONAL_WF,
    REPRO_ALT_GEN_P1,
    REPRO_ALT_GEN_P0,
    OLD_SURV_WF,
)


@dataclass
class HaploidWF(ReproductionBase):
    """Reproduction mode based on Wright-Fisher model with clonal
     haploid individuals."""
    
    pop_size: int #number of haploid individuals
    _gens_per_lifecycle: int = 1
    sexes: bool = False  # not yet used?

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """

        """Fill self.model.map with SLiM script snippets."""

        # methods inherited from parent Reproduction Base
        self._write_trees_file()

        self._define_subpopulations()
        self._add_initialize_constants()
        self._add_mode_scripts()
        self._add_survival_script()
        

    def _define_subpopulations(self):
        """Add a single diploid population. See NonWrightFisher for comparison."""
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts="")
        else:
            self.model.first(
                time=1,
                scripts="sim.addSubpop('p1', K, haploid=T);",
                comment="define starting haploid population.",
            )

    def _add_mode_scripts(self):
        """fitness and mating of diploid population."""

        self.model.repro(
            population="p1",
            scripts= REPRO_HAPLOID_WF,
            comment="haploid random mating; mating success weighted by fitness."
        )

    def _add_initialize_constants(self):
        """Add defineConstant calls to init for new variables."""
        metadata_dict = {
            'model': "shadie haploid WF with sex",
            'length': self.model.sim_time,
            'spo_pop_size': "NA",
            'gam_pop_size': self.pop_size,
            'spo_mutation_rate': self.model.metadata['mutation_rate'],
            'recombination_rate': self.model.metadata['recomb_rate']
        }

        self.model.map["initialize"][0]['constants']["K"] = self.pop_size
        self.model.map["initialize"][0]['simglobals']["METADATA"] = metadata_dict

    def _add_survival_script(self):
        """
        Defines the late() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        self.model.survival(
            population=None,
            scripts="return (individual.age == 0);",
            comment="non-overlapping generations",
        )

@dataclass
class ClonalHaploidWF(ReproductionBase):
    """Reproduction mode based on Wright-Fisher model with clonal
     haploid individuals."""
    
    pop_size: int #number of haploid individuals
    _gens_per_lifecycle: int = 1
    sexes: bool = False  # not yet used?

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        self._define_subpopulations()
        self._add_initialize_constants()
        self._add_scripts()
        self._add_survival_script()
        self._write_trees_file()

    def _define_subpopulations(self):
        """Add a single diploid population. See NonWrightFisher for comparison."""
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts="")
        else:
            self.model.first(
                time=1,
                scripts="sim.addSubpop('p1', K, haploid=T);",
                comment="define starting haploid population.",
            )

    def _add_scripts(self):
        """fitness and mating of diploid population."""

        self.model.repro(
            population="p1",
            scripts= REPRO_CLONAL_WF,
            comment="clonal reproduction; reproduction success weighted by fitness."
        )

    def _add_initialize_constants(self):
        """Add defineConstant calls to init for new variables."""
        metadata_dict = {
            'model': "shadie haploid clonal WF",
            'length': self.model.sim_time,
            'spo_pop_size': "NA",
            'gam_pop_size': self.pop_size,
            'spo_mutation_rate': self.model.metadata['mutation_rate'],
            'recombination_rate': self.model.metadata['recomb_rate']
        }

        self.model.map["initialize"][0]['constants']["K"] = self.pop_size
        self.model.map["initialize"][0]['simglobals']["METADATA"] = metadata_dict

    def _add_survival_script(self):
        """
        Defines the late() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        self.model.survival(
            population=None,
            scripts="return (individual.age == 0);",
            comment="non-overlapping generations",
        )

@dataclass
class AltGenWF(ReproductionBase):
    """Reproduction mode based on Wright-Fisher model with clonal
     haploid individuals."""
    
    spo_pop_size: int #number of diploid individuals
    gam_pop_size: int #number of haploid individuals
    _gens_per_lifecycle: int = 2
    sexes: bool = False  # not yet used?

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
         # methods inherited from parent WrightFisher:
        self._define_subpopulations()
        self._add_initialize_constants()
        self._add_scripts()
        self._write_trees_file()

        #specific to Moran (remove survival; Moran had non-overlapping gens):
        self._add_survival_script()
        self._define_subpopulations()
        self._add_initialize_constants()
        self._add_scripts()
        self._write_trees_file()

    def _define_subpopulations(self):
        """Add a single diploid population. See NonWrightFisher for comparison."""
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts="")
        else:
            self.model.first(
                time=1,
                scripts=("sim.addSubpop('p1', SPO_POP_SIZE);\n"
                		 "sim.addSubpop('p0', 0);"),
                comment="define starting haploid population.",
            )

    def _add_scripts(self):
        """fitness and mating of diploid population."""

        self.model.repro(
            population="p1",
            scripts= REPRO_ALT_GEN_P1,
            comment="clonal random mating; mating success weighted by fitness."
        )

        self.model.repro(
            population="p0",
            scripts= REPRO_ALT_GEN_P0,
            comment="clonal random mating; mating success weighted by fitness."
        )

    def _add_initialize_constants(self):
        """Add defineConstant calls to init for new variables."""
        metadata_dict = {
            'model': "shadie haploid clonal WF",
            'length': self.model.sim_time,
            'spo_pop_size': self.spo_pop_size,
            'gam_pop_size': self.gam_pop_size,
            'spo_mutation_rate': self.model.metadata['mutation_rate'],
            'recombination_rate': self.model.metadata['recomb_rate']
        }

        self.model.map["initialize"][0]['constants']["SPO_POP_SIZE"] = self.spo_pop_size
        self.model.map["initialize"][0]['constants']["GAM_POP_SIZE"] = self.gam_pop_size
        self.model.map["initialize"][0]['simglobals']["METADATA"] = metadata_dict

    def _add_survival_script(self):
        """
        Defines the late() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        self.model.survival(
            population=None,
            scripts="return (individual.age == 0);",
            comment="non-overlapping generations",
        )

    def _add_muteffect_script(self):
        idx = 6
        # iterate over MutationTypes
        for mut in self.model.chromosome.mutations:
            if mut._expr != "None":

                # refer to mutations by s{idx}
                idx += 1
                sidx = str("s" + str(idx))

                # add mutEffect callback function (e.g., s5 mutEffect(m1) {...})
                # for each MutationType. This callback will be activated or
                # deactivated (below) by early scripts based on whether
                # it is the haploid or diploid subpopulation's generation.
                if mut._expr == "haploid":
                    self.model.muteffect(
                        idx = None,
                        mutation = mut.name,
                        scripts = HAP_MUT_FITNESS,
                        comment = "mutation only expressed in haploid"
                        )
                elif mut._expr == "diploid":
                    self.model.muteffect(
                        idx = None,
                        mutation = mut.name,
                        scripts = DIP_MUT_FITNESS,
                        comment = "mutation only expressed in diploid"
                        )
                elif mut._expr == "None":
                    pass
                else:
                    print("Differental expression must be set to 'haploid'"
                        "or 'diploid")


@dataclass
class Moran(WrightFisher):
    """Reproduction mode based on Wright-Fisher model."""
    pop_size: int
    sexes: bool = False  # not yet used?
    _gens_per_lifecycle: int = 1

    def run(self):
        """
        Updates self.model.map with new component scripts for running
        life history and reproduction based on input args.
        """
        # methods inherited from parent WrightFisher:
        self._define_subpopulations()
        self._add_initialize_constants()
        self._add_scripts()
        self._write_trees_file()

        #specific to Moran (remove survival; Moran had non-overlapping gens):
        self._add_survival_script()

    def _add_survival_script(self):
        pass


if __name__ == "__main__":

    import shadie

    # define mutation types
    m0 = shadie.mtype(0.5, 'n', [0, 0.4])
    m1 = shadie.mtype(0.5, 'g', [0.8, 0.75])
    m2 = shadie.mtype(0.5, 'g', [0.8, 0.75], affects_haploid = False)
    m3 = shadie.mtype(0.5, 'n', [0, 0.4], affects_diploid = False)

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
        mod.reproduction.wright_fisher_haploid_clonal(
            pop_size = 500,
            #spo_pop_size=100,
            #gam_pop_size=100,
        )
    print(mod.script)
    #mod.write("/tmp/slim.slim")
    #mod.run(seed=123)
