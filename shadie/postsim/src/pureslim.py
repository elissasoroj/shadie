#!/usr/bin/env python

"""
A returned object class from a shadie simulation call.
"""

from typing import List, Optional, Iterable, Union
import numpy as np
import pyslim
import tskit
import msprime
import pandas as pd
import scipy.stats
from loguru import logger

import toyplot
from toytree.utils import toytree_sequence, ScrollableCanvas
from shadie.chromosome.src.classes import ChromosomeBase

class PureSlim:
    """"Loads two TreeSequence files with ancestral burn-in
    from a pure SLiM simulation
    """
    def __init__(
        self,
        tree_file: str,
        seed: Optional[int]=None,
        **kwargs,
        ):

        # hidden attributes
        self._tree_sequence = tskit.load(tree_file)

        # attributes to be parsed from the slim metadata
        self.generations: int=0
        """The SLiM simulated length of time in diploid generations AFTER
        the ancestral burn in (branch gens only)"""
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
        self._remove_null_population_and_nodes()

    def _extract_metadata(self):
        """Extract self attributes from shadie .trees file metadata.

        TODO: can more of this be saved in SLiM metadata?
        """
        gens = self._tree_sequence.metadata["SLiM"]["generation"]
        self.generations = gens
        assert self.popsize, "popsize not found in metadata; must enter a popsize arg."
        assert self.mut, "mut not found in metadata; must enter a mut arg."
        assert self.recomb, "recomb not found in metadata; must enter a recomb arg."

    def _remove_null_population_and_nodes(self):
        """Call tskit simplify function to remove null pop.

        There is a null population in shadie simulations because we
        define an alternation of generations with two alternating
        subpopulations. At the final generation of shadie SLiMulation
        the generation is even, and so we ...
        """
        # set population=0 for all nodes in each ts. Nodes from the
        # diploid sub-generation are currently labeled as population=1.
        

        # tskit tables are immutable, but we can modify a copy of
        # the table and use load_tables to make a new ts from it.
        tables = self._tree_sequence.tables

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
        self.tree_sequence = mod_tree_seq.simplify(
            samples=nodes_in_edge_table,
            keep_input_roots=True,
            keep_unary_in_individuals=True
        )

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
            sample_1 = samples[tts.sample[0]:]

            stats = pd.Series(
                index=["theta_0", "theta_1", "Fst_01", "Dist_01", "D_Taj_0", "D_Taj_1"],
                name=str(rep),
                data=[
                    tts.tree_sequence.diversity(sample_0),
                    tts.tree_sequence.diversity(sample_1),
                    tts.tree_sequence.Fst([sample_0, sample_1]),
                    tts.tree_sequence.divergence([sample_0, sample_1]),
                    tts.tree_sequence.Tajimas_D(sample_0),
                    tts.tree_sequence.Tajimas_D(sample_1),
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
        if self.chromosome is not None:
            self.chromosome.draw(axes=ax0)
            # ax0.y.show = False
            # ax0.x.ticks.show = True
            ax0.x.ticks.near = 5
            ax0.x.ticks.far = 0

        # add generation line showing where SLiM simulation ended.
        # only add if not much higher than highest tree height.
        # thsi could be faster by not building all these trees...
        top_root = max([i.treenode.height for i in tts][start:start+max_trees])
        if show_generation_line:
            if self.generations < top_root + top_root * 0.1:
                axes.hlines(
                    self.generations,
                    style={"stroke-dasharray": "4,2", "stroke-opacity": 0.4}
                )

        return canvas, axes, mark

    def draw_stats(
        self,
        stat: str="diversity",
        window_size: int=20_000,
        sample: Union[int, Iterable[int]]=6,
        reps: int=1,
        seed: Optional[int]=None,
        ):
        """Return a toyplot drawing of a statistic across the genome.
        
        If reps > 1 the measurement is repeated on multiple sets of 
        random samples of size `sample`, and the returned statistic
        is the mean with +/- 1 stdev shown. 

        Parameters
        ----------
        """
        # select a supported statistic to measure
        if stat == "diversity":
            func = self.tree_sequence.diversity
        else:
            raise NotImplementedError(f"stat {stat} on the TODO list...")

        # repeat measurement over many random sampled replicates
        rng = np.random.default_rng(seed)
        rep_values = []
        for _ in range(reps):
            samples = rng.choice(self.tree_sequence.samples(), sample, replace=False)
            values = func(
                sample_sets=samples,
                windows=np.linspace(
                    start=0, 
                    stop=self.tree_sequence.sequence_length, 
                    num=round(self.tree_sequence.sequence_length / window_size)
                )
            )
            rep_values.append(values)
        
        # get mean and std
        means = np.array(rep_values).mean(axis=0)
        stds = np.array(rep_values).mean(axis=0)        

        # draw canvas...
        canvas, axes, mark  = toyplot.fill(
            means, height=300, width=500, opacity=0.5, margin=(60, 50, 50, 80)
        )

        # style axes
        axes.x.ticks.show = True
        axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
        axes.y.ticks.labels.angle = -90
        axes.y.ticks.show = True
        axes.y.ticks.locator = toyplot.locator.Extended(only_inside=True, count=5)        
        axes.label.offset = 20
        axes.label.text = f"{stat} in {int(window_size / 1000)}kb windows"
        return canvas, axes, mark




class PureSlim_TwoPops:
    """"Loads and merges two TreeSequence files with ancestral burn-in
    from a pure SLiM simulation
    
    Two runs were initiated from an ancestral SLiM burn-in .trees file

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
        """The SLiM simulated length of time in diploid generations AFTER
        the ancestral burn in (branch gens only)"""
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
        self._union()

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

    def _union(self):
        """Calls tskit union.
        ...
        """
        self._remove_null_population_and_nodes()
        self._merge_ts_pops()
        self._report_ninds()

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
        if self._nts > 2:
            raise ValueError("you cannot enter >2 tree sequences.") 
        
        node_map= self._match_nodes(other=self._tree_sequences[0],
            ts=self._tree_sequences[1],
            split_time= int(1+(2*self.generations)))

        tsu = self._tree_sequences[1].union(self._tree_sequences[0],
            node_map, check_shared_equality=True)

        self.tree_sequence = pyslim.SlimTreeSequence(tsu)

    def _match_nodes(self, other, ts, split_time):
        """
        Given SLiM tree sequences `other` and `ts`, builds a numpy array with length
        `other.num_nodes` in which the indexes represent the node id in `other` and the
        entries represent the equivalent node id in `ts`. If a node in `other` has no
        equivalent in `ts`, then the entry takes the value `tskit.NULL`. The matching
        is done by comparing the IDs assigned by SLiM which are kept in the NodeTable
         metadata. Further, this matching of SLiM IDs is done for times (going 
         backward-in-time) greater than the specified `split_time`.
        """

        node_mapping = np.full(other.num_nodes, tskit.NULL)
        sids0 = np.array([n.metadata["slim_id"] for n in ts.nodes()])
        sids1 = np.array([n.metadata["slim_id"] for n in other.nodes()])
        alive_before_split1 = (other.tables.nodes.time >= split_time)
        is_1in0 = np.isin(sids1, sids0)
        both = np.logical_and(alive_before_split1, is_1in0)
        sorted_ids0 = np.argsort(sids0)
        matches = np.searchsorted(
            sids0,
            sids1[both],
            side='left',
            sorter=sorted_ids0)
        node_mapping[both] = sorted_ids0[matches]
        
        return node_mapping

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
    

    def alt_stats(self, sample:int = 10, seed = 123):
        inds_alive_in_pop0 = self.tree_sequence.individuals_alive_at(time=0, population=0)
        inds_alive_in_pop1 = self.tree_sequence.individuals_alive_at(time=0, population=1)

        rng = np.random.default_rng(seed)

        sample0 = rng.choice(inds_alive_in_pop0, sample, replace=False)
        sample1 = rng.choice(inds_alive_in_pop1, sample, replace=False)

        nodes0 = []
        nodes1 = []

        for i in sample0:
           nodes0.extend(self.tree_sequence.individual(i).nodes)
        for i in sample1:
           nodes1.extend(self.tree_sequence.individual(i).nodes)

        rng = np.random.default_rng(seed)
        data = []
        stats = pd.Series(
            index=["theta_0", "theta_1", "Fst_01", "Dist_01", "D_Taj_0", "D_Taj_1"],
            name=str(rep),
            data=[
                self.tree_sequence.diversity(nodes0),
                self.tree_sequence.diversity(nodes1),
                self.tree_sequence.Fst([nodes0, nodes1]),
                self.tree_sequence.divergence([nodes0, nodes1]),
                self.tree_sequence.Tajimas_D(nodes1),
                self.tree_sequence.Tajimas_D(nodes0),
            ],
            dtype=float,
        )
        data.append(stats)

        # concat to a dataframe
        data = pd.concat(data, axis=1).T
        
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
            sample_1 = samples[tts.sample[0]:]

            stats = pd.Series(
                index=["theta_0", "theta_1", "Fst_01", "Dist_01", "D_Taj_0", "D_Taj_1"],
                name=str(rep),
                data=[
                    tts.tree_sequence.diversity(sample_0),
                    tts.tree_sequence.diversity(sample_1),
                    tts.tree_sequence.Fst([sample_0, sample_1]),
                    tts.tree_sequence.divergence([sample_0, sample_1]),
                    tts.tree_sequence.Tajimas_D(sample_0),
                    tts.tree_sequence.Tajimas_D(sample_1),
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
        if self.chromosome is not None:
            self.chromosome.draw(axes=ax0)
            # ax0.y.show = False
            # ax0.x.ticks.show = True
            ax0.x.ticks.near = 5
            ax0.x.ticks.far = 0

        # add generation line showing where SLiM simulation ended.
        # only add if not much higher than highest tree height.
        # thsi could be faster by not building all these trees...
        top_root = max([i.treenode.height for i in tts][start:start+max_trees])
        if show_generation_line:
            if self.generations < top_root + top_root * 0.1:
                axes.hlines(
                    self.generations,
                    style={"stroke-dasharray": "4,2", "stroke-opacity": 0.4}
                )

        return canvas, axes, mark

    def draw_stats(
        self,
        stat: str="diversity",
        window_size: int=20_000,
        sample: Union[int, Iterable[int]]=6,
        reps: int=1,
        seed: Optional[int]=None,
        ):
        """Return a toyplot drawing of a statistic across the genome.
        
        If reps > 1 the measurement is repeated on multiple sets of 
        random samples of size `sample`, and the returned statistic
        is the mean with +/- 1 stdev shown. 

        Parameters
        ----------
        """
        # select a supported statistic to measure
        if stat == "diversity":
            func = self.tree_sequence.diversity
        else:
            raise NotImplementedError(f"stat {stat} on the TODO list...")

        # repeat measurement over many random sampled replicates
        rng = np.random.default_rng(seed)
        rep_values = []
        for _ in range(reps):
            samples = rng.choice(self.tree_sequence.samples(), sample, replace=False)
            values = func(
                sample_sets=samples,
                windows=np.linspace(
                    start=0, 
                    stop=self.tree_sequence.sequence_length, 
                    num=round(self.tree_sequence.sequence_length / window_size)
                )
            )
            rep_values.append(values)
        
        # get mean and std
        means = np.array(rep_values).mean(axis=0)
        stds = np.array(rep_values).mean(axis=0)        

        # draw canvas...
        canvas, axes, mark  = toyplot.fill(
            means, height=300, width=500, opacity=0.5, margin=(60, 50, 50, 80)
        )

        # style axes
        axes.x.ticks.show = True
        axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
        axes.y.ticks.labels.angle = -90
        axes.y.ticks.show = True
        axes.y.ticks.locator = toyplot.locator.Extended(only_inside=True, count=5)        
        axes.label.offset = 20
        axes.label.text = f"{stat} in {int(window_size / 1000)}kb windows"
        return canvas, axes, mark
