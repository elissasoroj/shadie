#!/usr/bin/env python

"""
A returned object class from a shadie simulation call.
"""

from typing import List, Optional, Union, Iterable
import numpy as np
import pyslim
import tskit
import msprime
import pandas as pd
import scipy.stats
from loguru import logger

from toytree.utils.tree_sequence.src.toytree_sequence import ToyTreeSequence
from toytree.utils.utils import ScrollableCanvas
from shadie.chromosome.src.classes import ChromosomeBase


class TwoSims:
    """"Loads and merges two TreeSequence files from SLiM.

    SliM has been run for two populations with the same genome type
    for T generations starting from different seeds and perhaps
    with different selection or demographic histories. Each ts
    includes N haploid individuals at time_start and time_end.

    SLiM was run with the argument keep_input_roots=True so that it
    preserves the first-generation ancestors so they can be used for
    recapitation.
    """
    def __init__(
        self,
        tree_files: List[str],
        seed: Optional[int]=None,
        **kwargs,
        ):

        # hidden attributes
        self._tree_files: List[str] = tree_files
        self._tree_sequences = [pyslim.load(i) for i in self._tree_files]
        self._nts: int = len(self._tree_files)

        # attributes to be parsed from the slim metadata
        self.generations: int=0
        """The SLiM simulated length of time in diploid generations."""
        self.popsize: int=kwargs.get("popsize")
        """The SLiM simulated diploid carrying capacity"""
        self.recomb: float=kwargs.get("recomb")
        """The recombination rate to use for recapitation."""
        self.mut: float=kwargs.get("mut")
        """The mutation rate to use for recapitated treesequence."""
        self.chromosome: ChromosomeBase=kwargs.get("chromosome")
        """The shadie.Chromosome class representing the SLiM genome."""
        self.rng: np.random.Generator=np.random.default_rng(seed)

        # new attributes built as results
        self.tree_sequence: pyslim.SlimTreeSequence=None
        """A SlimTreeSequence that has been recapitated and mutated."""

        # try to fill attributes by extracting metadata from the tree files.
        self._extract_metadata()
        self._tskit_complete()

    def _extract_metadata(self):
        """Extract self attributes from shadie .trees file metadata.

        TODO: can more of this be saved in SLiM metadata?
        """
        gens = [i.metadata["SLiM"]["generation"] for i in self._tree_sequences]
        assert len(set(gens)) == 1, ("simulations must be same length (gens).")
        self.generations = gens[0]
        assert self.popsize, "popsize not found in metadata; must enter a popsize arg."
        assert self.mut, "mut not found in metadata; must enter a mut arg."
        assert self.recomb, "recomb not found in metadata; must enter a recomb arg."

    def _tskit_complete(self):
        """Calls tskit union, simplify, mutate to complete neutral sim.
        ...
        """
        self._remove_null_population_and_nodes()
        self._merge_ts_pops()
        self._report_ninds()
        self._divide_tree_height()
        self._recapitate()
        self._mutate()

    def _remove_null_population_and_nodes(self):
        """Call tskit simplify function to remove null pop.

        There is a null population in shadie simulations because we
        define an alternation of generations with two alternating
        subpopulations. At the final generation of shadie SLiMulation
        the generation is even, and so we ...
        """
        # set population=0 for all nodes in each ts. Nodes from the
        # diploid sub-generation are currently labeled as population=1.
        for idx, treeseq in enumerate(self._tree_sequences):

            # tskit tables are immutable, but we can modify a copy of
            # the table and use load_tables to make a new ts from it.
            tables = treeseq.tables

            # modify the tables to set population to 0 for all
            nnodes = tables.nodes.num_rows
            tables.nodes.population = np.zeros(nnodes, dtype=np.int32)

            # drop nodes that are not connected to anything. This includes
            # the pseudo-nodes representing half of the haploid populations.
            nodes_in_edge_table = list(
                set(tables.edges.parent).union(tables.edges.child)
            )

            # reload treeseq FROM modified tables
            mod_tree_seq = pyslim.load_tables(tables)

            # remove the empty population (p1) by using simplify, which will
            # find that there are no longer any nodes in population=1. This
            # does not remove any Nodes, but it does remove a population.
            # https://tskit.dev/tskit/docs/stable/_modules/tskit/tables.html
            self._tree_sequences[idx] = mod_tree_seq.simplify(
                samples=nodes_in_edge_table,
                keep_input_roots=True,
                keep_unary_in_individuals=True
            )

    def _merge_ts_pops(self):
        """Merge two separate sims into a single ts with 2 pops.

        Merges with union and re-loads the ts as a SlimTreeSequence.
        Adds the non-shared portions of ts1 to ts0. Since they have
        no shared portion, we enter NULL for the `node mapping`, and
        add_population True sets new nodes to a new population.
        """
        # check the number of simulations for what to do.
        if self._nts == 1:
            self.tree_sequence = pyslim.SlimTreeSequence(self._tree_sequences[0])
            return
        if self._nts > 2:
            raise ValueError("you cannot enter >2 tree sequences.") 
        # Merge two tree sequences
        ts0 = self._tree_sequences[0]
        ts1 = self._tree_sequences[1]
        merged_ts = ts0.union(
            ts1,
            node_mapping=[tskit.NULL for i in range(ts1.num_nodes)],
            add_populations=True,
        )
        self.tree_sequence = pyslim.SlimTreeSequence(merged_ts)

    def _report_ninds(self):
        """Report number of inds in each population."""
        treeseq = self.tree_sequence
        # get array of ids of individuals (sets of 2 nodes) alive at
        # 0 generations ago (diploids) [time=1 also has nodes].
        # e.g., [0, 1, 2, ... 2000, 2001, 2002]
        inds_alive_in_pop0 = treeseq.individuals_alive_at(time=0, population=0)
        inds_alive_in_pop1 = treeseq.individuals_alive_at(time=0, population=1)
        npop0 = len(inds_alive_in_pop0)
        npop1 = len(inds_alive_in_pop1)
        logger.info(f"inds alive at time=0; simpop0={npop0}, simpop1={npop1}")

    def _recapitate(self):
        """Merge pops backwards in time and simulate ancestry.
        """
        # recapitate: ts is passed to sim_ancestry as 'initial_state'.
        # this automatically merges everyone into new ancestral pop.
        self.tree_sequence = pyslim.recapitate(
            ts = self.tree_sequence,
            recombination_rate = self.recomb,
            ancestral_Ne=self.popsize,
            random_seed=self.rng.integers(2**31),
        )

    def _divide_tree_height(self):
        """Divide all time measurements by half in ts tables.

        The number of generations in a shadie simulation is 2X that of
        a normal simulations since a single generation represents just
        the haploid or diploid phase of a normal generation. To make
        the trees normal again we divide time by 2X.

        Mutation and recombination rates are handled separately.
        """
        # set metadata
        self.tree_sequence.metadata["SLiM"]["generation"] = int(
            self.tree_sequence.metadata["SLiM"]["generation"] / 2)

        # tskit tables are immutable, get a copy
        tables = self.tree_sequence.tables
        tables.mutations.time = tables.mutations.time / 2.
        tables.nodes.time = tables.nodes.time / 2.

        # reload treeseq from modified tables
        self.tree_sequence = pyslim.load_tables(tables)

    def _mutate(self):
        """Mutatates the recapitated TreeSequence.

        This applies a mutation model to edges of the tree sequence.
        Does it know which regions to mutate or not mutate? For example,
        all recapitated edges should be mutated, but also the neutral
        genomic regions of the SLiM time frame should be mutated.
        """
        # logger report before adding mutations
        self._report_mutations(allow_m0=False)

        # add mutations
        self.tree_sequence = msprime.sim_mutations(
            self.tree_sequence,
            rate=self.mut,
            random_seed=self.rng.integers(2**31),
            keep=True,  # whether to keep existing mutations.
            model=msprime.SLiMMutationModel(type=0),
        )
        self.tree_sequence = pyslim.SlimTreeSequence(self.tree_sequence)

        # logger report after adding mutations
        self._report_mutations(allow_m0=True)

    def _report_mutations(self, allow_m0: bool=True):
        """Report to logger nmutations, mtypes, and check for m0 type."""
        # check type m0 is not used
        mut_types = set(
            mutlist['mutation_type']
            for mut in self.tree_sequence.mutations()
            for mutlist in mut.metadata['mutation_list']
        )
        if not allow_m0:
            assert 0 not in mut_types, "m0 mutation types already present."

        # report to logger the existing mutations
        logger.info(
            f"Keeping {self.tree_sequence.num_mutations} existing "
            f"mutations of type(s) {mut_types}.")

    def stats(
        self,
        sample: Union[int, Iterable[int]]=10,
        seed: Optional[int]=None,
        reps: int=10
        ):
        """Calculate statistics summary on mutated TreeSequence.
        
        Returns a dataframe with several statistics calculated and
        summarized from replicate random sampling.

        Parameters
        ----------
        sample: int or Iterable of ints
            The number of tips to randomly sample from each population.
        seed: int
            A seed for random sampling.
        reps: int
            Number of replicate times to random sample tips and 
            calculate statistics.

        Returns
        -------
        pandas.DataFrame
            A dataframe with mean and 95% confidence intervals.
        """
        rng = np.random.default_rng(seed)
        data = []

        # get a list of Series
        for rep in range(reps):
            seed = rng.integers(2**31)
            tts = ToyTreeSequence(self.tree_sequence, sample=sample, seed=seed)
            samples = np.arange(tts.sample[0] + tts.sample[1])
            sample_0 = samples[:tts.sample[0]]
            sample_0_nodes = []
            for i in sample_0:
               sample_0_nodes.extend(self.tree_sequence.individual(i).nodes)
            
            sample_1 = samples[tts.sample[0]:]
            sample_1_nodes = []
            for i in sample_1:
               sample_1_nodes.extend(self.tree_sequence.individual(i).nodes)

            stats = pd.Series(
                index=["theta_0", "theta_1", "Fst_01", "Dist_01", "D_Taj_0", "D_Taj_1"],
                name=str(rep),
                data=[
                    tts.tree_sequence.diversity(sample_0_nodes),
                    tts.tree_sequence.diversity(sample_1_nodes),
                    tts.tree_sequence.Fst([sample_0_nodes, sample_1_nodes]),
                    tts.tree_sequence.divergence([sample_0_nodes, sample_1_nodes]),
                    tts.tree_sequence.Tajimas_D(sample_0_nodes),
                    tts.tree_sequence.Tajimas_D(sample_1_nodes),
                ],
                dtype=float,
            )
            data.append(stats)

        # concat to a dataframe
        data = pd.concat(data, axis=1).T

        # get 95% confidence intervals
        confs = []
        for stat in data.columns:
            mean_val = np.mean(data[stat])
            low, high = scipy.stats.t.interval(
                alpha=0.95,
                df=len(data[stat]) - 1,
                loc=mean_val,
                scale=scipy.stats.sem(data[stat]),
            )
            confs.append((mean_val, low, high))

        # reshape into a dataframe
        data = pd.DataFrame(
            columns=["mean", "CI_5%", "CI_95%"],
            index=data.columns,
            data=np.vstack(confs),
        )
        return data


    def draw_tree(
        self,
        idx: int=None,
        site: int=None,
        sample: Union[int, Iterable[int]]=10,
        seed=None,
        show_label=True,
        show_generation_line=True,
        **kwargs):
        """Returns a toytree drawing for a random sample of tips.

        The tree drawing will include mutations as marks on edges
        with mutationTypes colored differently. The tips are re-labeled
        indicating {population-id}-{random-sample-id}.

        Parameters
        ----------
        ...
        """
        # load as a ToyTreeSequence
        tts = ToyTreeSequence(self.tree_sequence, sample=sample, seed=seed)

        # draw tree and mutations with a pre-set style
        base_style = {
            'scale_bar': True,
            'width': 400,
            'height': 400,
        }
        base_style.update(kwargs)
        canvas, axes, marks = tts.draw_tree(
            idx=idx, site=site, show_label=show_label, **base_style)

        # add fill to show the SLiM portion of the simulation.
        if show_generation_line:
            style = {"stroke-dasharray": "4,2", "stroke-width": 2, "stroke": "grey"}
            if marks[0].layout == "u":
                axes.hlines(-self.generations, style=style)
            elif marks[0].layout == "d":
                axes.hlines(self.generations, style=style)
            elif marks[0].layout == "r":
                axes.vlines(-self.generations, style=style)
            elif marks[0].layout == "l":
                axes.vlines(self.generations, style=style)
        return canvas, axes, marks


    def draw_tree_sequence(
        self,
        start: int=0,
        max_trees: int=10,
        sample: Union[int,Iterable[int]]=10,
        seed: Optional[int]=None,
        height: Optional[int]=None,
        width: Optional[int]=None,
        show_generation_line: bool=True,
        scrollable: bool=True,
        axes: Optional['toyplot.coordinates.Cartesian']=None,
        **kwargs,
        ):
        """Return a ToyTreeSequence drawing.

        """
        # load as a ToyTreeSequence
        tts = ToyTreeSequence(self.tree_sequence, sample=sample, seed=seed)

        # get auto-dimensions from tree size
        height = height if height is not None else 325
        height += 100
        width = max(300, (
            width if width is not None else 
            15 * tts.at_index(0).ntips * min(max_trees, len(tts))
        ))        

        # create an axis for the chromosome. +100 height for chrom.
        canvas = ScrollableCanvas(height=height, width=width)
        ax0 = canvas.cartesian(bounds=(50, -50, 50, 75), padding=5)
        ax1 = canvas.cartesian(bounds=(50, -50, 100, -50))        

        # draw tree sequence
        canvas, axes, mark = tts.draw_tree_sequence(
            start=start,
            max_trees=max_trees,
            axes=ax1,
            **kwargs,
        )

        # add chromosome to top axis
        self.chromosome.draw(axes=ax0)
        # ax0.y.show = False
        # ax0.x.ticks.show = True
        ax0.x.ticks.near = 5
        ax0.x.ticks.far = 0
        # ax0.fill(
            # [0, 1], [0, 0], [1, 1], 
            # color='green',
        # )
        # ax0.x.ticks.locator = toyplot.locator.Extended(count=8)


        # ntrees = len(tts)
        # tmax = start + min(ntrees, max_trees)
        # breaks = tts.tree_sequence.breakpoints(True)[start: tmax + 1]
        # cend = breaks[-1]
        # print(cend)

        # cdat = self.chromosome.data.loc[start:tmax]



        # add generation line showing where SLiM simulation ended.
        if show_generation_line:
            axes.hlines(
                self.generations,
                style={"stroke-dasharray": "4,2", "stroke-opacity": 0.4}
            )

        return canvas, axes, mark