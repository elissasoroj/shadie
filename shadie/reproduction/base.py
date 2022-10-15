#!/usr/bin/env python

"""
Base classes for reproduction.

ReproductionBase -> NonWrightFisher -> BryophyteBase
                    WrightFisher       PteridophyteBase
                                       etc.
"""

from dataclasses import dataclass
import pyslim
from shadie.reproduction.scripts import (
    SUBSTITUTION,
    SUB_MUTS,
    P0_FITNESS_SCALE_DEFAULT,
    P1_FITNESS_SCALE_DEFAULT,
    EARLY,
    EARLY_WITH_GAM_K,
    WF_REPRO,
    HAP_MUT_FITNESS,
    DIP_MUT_FITNESS
)


@dataclass
class ReproductionBase:
    """Reproduction code block generation BaseClass.

    All Reproduction subclasses will store the Model, lineage,
    mode, and a set of scripts stored as a dictionary. Users will
    only interact with subclasses of ReproductionBase, generated by
    using the ReproductionApi accessible from the Model object.

    The NonWFBase is a superclass used for organism models. The
    WFBase superclass is used mainly for testing or comparison to
    standard WF model simulations.

    Example
    -------
    with shadie.Model() as model:
        model.initialize(...)
        model.early(...)
        model.reproduction.bryophyte(...)
    """
    model: 'shadie.Model'

    def _write_trees_file(self):
        """adds late() call to save and write .trees file.

        All shadie reproduction classes write a .trees file in a late()
        call, but the time at which to write it varies depending on
        whether the start point was loaded from a previous file.
        """
        # get time AFTER the last even generation.
        endtime = int(self.model.sim_time + 1)

        # calculate end based on this sim AND the loaded parent sim.
        if self.model.metadata['file_in']:
            ts_start = tskit.load(self.model.metadata['file_in'])
            sim_start = ts_start.max_root_time
            resched_end = int(endtime + sim_start)
            self.model.late(
                time=resched_end,
                scripts=[
                    "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)",
                    f"sim.treeSeqOutput('{self.model.metadata['file_out']}', metadata = METADATA)"],
                comment="end of sim; save .trees file",
            )
        # write output at last generation of this simulation.
        else:
            self.model.late(
                time=endtime,
                scripts=[
                    "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)",
                    f"sim.treeSeqOutput('{self.model.metadata['file_out']}', metadata = METADATA)"],
                comment="end of sim; save .trees file",
            )


@dataclass
class NonWrightFisher(ReproductionBase):
    """Reproduction mode based on NON Wright-Fisher model.

    This is a subclass that is extended for organism specific
    reproduction models in shadie.reproduction. All NonWF models
    include alternation of generations (p0 and p1 subpops). The
    alternative is to implement a WF model.
    """
    def _set_gametophyte_k(self):
        """Sets a carrying capacity for gametophyte holding pop (during p1
        generation, to avoid lagging in the simulation. Automatically sets
        to 10x user-defined popsize
        """
        if not self.gam_k:
            self.gam_k = 10*self.gam_pop_size

    def _define_subpopulations(self):
        """add haploid and diploid life stages as subpopulations."""
        if self.model.metadata['file_in']:
            self.model._read_from_file(tag_scripts =["p1.individuals.tag=3;", 
                "tags = rbinom(1, p0.individualCount, 0.5);", "p0.individuals.tag = tags;"])
        else:
            self.model.early(
                time=1,
                scripts=[
                    "sim.addSubpop('p1', SPO_POP_SIZE)",
                    "sim.addSubpop('p0', 0)",
                    "p1.individuals.tag = 3",],
                comment="define subpops: p1=diploid sporophytes, p0=haploid gametophytes",
            )

    def _add_initialize_globals(self):
        """Add defineGlobal calls to init variables.

        When this is called by different superclasses that have
        different attributes unique to each it stores only their
        unique set of attributes into a metadata dictionary that
        will be saved into the .trees file at the end of the sim.
        Includes parent attrs like model.
        """
        # exclude parent class attributes
        exclude = ["_substitution_str", "model",
                    "_p0activate_str", "_p0deactivate_str",
                    "_p1activate_str", "_p1deactivate_str"]
        asdict = {
            i: j for (i, j) in self.__dict__.items()
            if i not in exclude
        }
        self.model.map["initialize"][0]['simglobals']['METADATA'] = asdict

    def _add_initialize_constants(self):
        """Add defineConstant calls to init variables.

        When this is called by different superclasses that have
        different attributes unique to each it stores only their
        unique set of attributes. Excludes parent attrs like model.
        """
        # exclude parent class attributes
        exclude = ["lineage", "mode", "model", "_substitution_str", 
                    "_p0activate_str", "_p0deactivate_str",
                    "_p1activate_str", "_p1deactivate_str"]
        asdict = {
            i: j for (i, j) in self.__dict__.items()
            if i not in exclude
        }
        self.model.map["initialize"][0]['constants'].update(asdict)

    def _add_alternation_of_generations(self):
        """Alternation of generations scripts.

        This writes fitness, early, and late functions that activate
        or deactivate fitness effects of mutations in alternating
        generations.
        """
        # add fitness callback for gametophytes based on MutationTypes
        # in the model.chromosome.
        # this will map to sx-sy survival callbacks.
        idx = 6
        p0activate_scripts = []
        p0deactivate_scripts = []
        p1activate_scripts = []
        p1deactivate_scripts = []
        substitutions = []

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

                # add reference to this mutation to be added to a late call
                # for checking whether a mutation has become a substitution.
                #CHECK - SHOULD NOT BE NECESSARY ANYMORE
                # sub_muts = SUB_MUTS.format(idx=sidx, mut=mut.name).lstrip()
                # substitutions.append(sub_muts)

        # insert references to fitness callbacks into an early script
        # that will alternately activate or deactivate them on
        # alternating generations to only apply to gameto or sporo.
        p0activate_str = "\n        ".join(p0activate_scripts)
        p0deactivate_str = "\n        ".join(p0deactivate_scripts)
        p1activate_str = "\n        ".join(p1activate_scripts)
        p1deactivate_str = "\n        ".join(p1deactivate_scripts)

        #save activate and deactivate scripts for later
        self._p0activate_str = p0activate_str
        self._p0deactivate_str = p0deactivate_str
        self._p1activate_str = p1activate_str
        self._p1deactivate_str = p1deactivate_str

        #CHECK - SHOULD NOT BE NECESSARY ANYMORE
        # # insert the substitution-checking scripts into larger context
        # substitution_str = "\n    ".join(substitutions)
        # #save subsitutions for late caldl in model-specific scripts
        # self._substitution_str = substitution_str

    def _add_early_script(self):
        """
        Defines the early() callbacks for each gen.
        This overrides the NonWrightFisher class function of same name.
        """
        early_script = (EARLY_WITH_GAM_K.format(
            p0_fitnessScaling= P0_FITNESS_SCALE_DEFAULT,
            p1_fitnessScaling= P1_FITNESS_SCALE_DEFAULT,
            p0activate= self._p0activate_str,
            p0deactivate= self._p0deactivate_str,
            p1deactivate= self._p1deactivate_str,
            p1activate= self._p1activate_str,
            )
        )

        self.model.early(
            time=None,
            scripts=early_script,
            comment="alternation of generations",
        )


@dataclass
class WrightFisher(ReproductionBase):
    """Reproduction mode based on Wright-Fisher model."""
    pop_size: int
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
            self.model.early(
                time=1,
                scripts="sim.addSubpop('p1', K);",
                comment="define starting diploid population.",
            )

    def _add_scripts(self):
        """fitness and mating of diploid population."""

        self.model.repro(
            population="p1",
            scripts= WF_REPRO,
            comment="hermaphroditic random mating."
        )

    def _add_initialize_constants(self):
        """Add defineConstant calls to init for new variables."""
        metadata_dict = {
            'model': "shadie WF",
            'length': self.model.sim_time,
            'spo_pop_size': self.pop_size,
            'gam_pop_size': "NA",
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


if __name__ == "__main__":

    import shadie

    # define mutation types
    m0 = shadie.mtype(0.5, 'n', 0, 0.4)
    m1 = shadie.mtype(0.5, 'g', 0.8, 0.75, diffexpr="diploid")

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

    print(m1._expr)

    with shadie.Model() as mod:
        mod.initialize(chromosome=chrom, sim_time=1000, #file_in = "/tmp/test.trees"
            )
        mod.reproduction.wright_fisher(pop_size=1000)
    print(mod.script)
    #mod.write("/tmp/slim.slim")
    #mod.run(binary="/usr/local/bin/slim")
