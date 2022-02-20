#!/usr/bin/env python

"""
A returned object class from a shadie simulation call.
"""

from typing import Optional, Union, Iterable, List
from dataclasses import dataclass, field
import numpy as np
import pyslim
# import tskit
import msprime
import toyplot
import pandas as pd
import tskit
import scipy.stats
from loguru import logger

from toytree.utils.src import toytree_sequence, ScrollableCanvas
from shadie.chromosome.src.classes import ChromosomeBase

logger = logger.bind(name='shadie')


class OneSim:
    """"Load and merge one TreeSequence files from SLiM.

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
        trees_file: str, #'tskit.trees.TreeSequence',
        ancestral_Ne: int,
        mut: float,
        recomb: float,
        chromosome: 'shadie.Chromosome',
        seed: Optional[int]=None,
        recapitate: bool=False,
        add_neutral_mutations: bool=False,
        ):

        # hidden attributes
        self.tree_sequence = pyslim.load(trees_file)
        """A SlimTreeSequence that has been recapitated and mutated."""

        # attributes to be parsed from the slim metadata
        self.generations: int=0
        """The SLiM simulated length of time in diploid generations."""
        self.ancestral_Ne: int=ancestral_Ne
        """The SLiM simulated diploid carrying capacity"""
        self.recomb: float=recomb
        """The recombination rate to use for recapitation."""
        self.mut: float=mut
        """The mutation rate to use for recapitated treesequence."""
        self.chromosome: ChromosomeBase=chromosome
        """The shadie.Chromosome class representing the SLiM genome."""
        self.rng: np.random.Generator=np.random.default_rng(seed)

        # try to fill attributes by extracting metadata from the tree files.
        self._extract_metadata()
        self._update_tables()
        if recapitate:
            self._recapitate()
        if add_neutral_mutations:
            self._mutate()

    def _extract_metadata(self):
        """Extract self attributes from shadie .trees file metadata.

        TODO: can more of this be saved in SLiM metadata?
        """
        self.generations = self.tree_sequence.metadata["SLiM"]["generation"]
        assert self.ancestral_Ne, "ancestral_Ne not found in metadata; must enter an ancestral_Ne arg."
        assert self.mut, "mut not found in metadata; must enter a mut arg."
        assert self.recomb, "recomb not found in metadata; must enter a recomb arg."

    def _update_tables(self):
        """Remove extra psuedopopulation nodes."""
        # get mutable tskit.TableCollection
        tables = self.tree_sequence.dump_tables()
        nnodes = tables.nodes.time.size

        # there is a null SLiM population (0) that doesnt really exist
        # and the actual poulation (1). So we set all to 0.
        tables.nodes.population = np.zeros(nnodes, dtype=np.int32)
        #meta = tables.metadata
        #meta["SLiM"]["generation"] = int(meta["SLiM"]["generation"] / 2.)
        #tables.metadata = meta
        #tables.nodes.time /= 2
        #tables.mutations.time /= 2.

        # drop nodes that are not connected to anything. This includes
        # the pseudo-nodes representing half of the haploid populations.
        nodes_in_edge_table = list(
            set(tables.edges.parent).union(tables.edges.child)
        )

        # remove the empty population nodes by using simplify, which 
        tables.simplify(
            samples=nodes_in_edge_table,
            keep_input_roots=True,
            filter_individuals=True,
            filter_populations=True,
            filter_sites=False,
        )

        # turn it back into a treesequence
        self.tree_sequence = pyslim.load_tables(tables)

    def _recapitate(self):
        """Merge pops backwards in time and simulate ancestry.
        """
        # recapitate: ts is passed to sim_ancestry as 'initial_state'.
        # this automatically merges everyone into new ancestral pop.
        self.tree_sequence = pyslim.recapitate(
            ts=self.tree_sequence,
            ancestral_Ne=self.ancestral_Ne,
            random_seed=self.rng.integers(2**31),
            recombination_rate=self.recomb,
        )

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
            tts = toytree_sequence(self.tree_sequence, sample=sample, seed=seed)
            samples = np.arange(tts.sample[0])

            stats = pd.Series(
                index=["theta", "D_Taj"],
                name=str(rep),
                data=[
                    tts.tree_sequence.diversity(samples),
                    tts.tree_sequence.Tajimas_D(samples),
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
        show_generation_line=False,
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
        tts = toytree_sequence(self.tree_sequence, sample=sample, seed=seed)

        # draw tree and mutations with a pre-set style
        base_style = {
            'scale_bar': True,
            'width': 300,
            'height': 300,
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
        # load as a ToyTreeSequence (TODO: make faster)
        tts = toytree_sequence(self.tree_sequence, sample=sample, seed=seed)

        # get auto-dimensions from tree size
        height = height if height is not None else 325
        height += 100
        width = max(300, (
            width if width is not None else 
            15 * tts.at_index(0).ntips * min(max_trees, len(tts))
        ))        

        # create an axis for the chromosome. +100 height for chrom.
        if scrollable:
            canvas = ScrollableCanvas(height=height, width=width)
        else:
            canvas = toyplot.Canvas(height=height, width=width)
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
        color: Optional[str]="lightseagreen"
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
            ndt = self.tree_sequence.tables.nodes
            mask = (ndt.population == 0) & (ndt.time == 0) & (ndt.flags == 1)
            arr = np.arange(mask.shape[0])[mask]
            size = min(arr.size, sample)
            samples = rng.choice(arr, size=size, replace=False)

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
        style = {"fill":str(color)}

        canvas, axes, mark  = toyplot.fill(
            means, height=300, width=500, opacity=0.5, margin=(60, 50, 50, 80), 
            style=style,
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


@dataclass
class TwoSims:
    trees_files: List[str]
    ancestral_ne: int
    mut: float
    recomb: float
    chromosome: 'shadie.Chromosome'
    seed: Optional[int]=None
    recapitate: bool=False
    add_neutral_mutations: bool=False

    _tree_sequences: List['tskit.trees.TreeSequence'] = field(default_factory=list, init=False)
    tree_sequence: 'tskit.trees.TreeSequence' = field(default=None, init=False)

    def __post_init__(self):
        self._update_tables()
        self._extract_metadata()
        self._match_nodes_and_merge()
        if self.recapitate:
            self._recapitate()
        if self.add_neutral_mutations:
            pass

    def _extract_metadata(self):
        """Extract self attributes from shadie .trees file metadata.
        TODO: can more of this be saved in SLiM metadata?
        """
        gens = [i.metadata["SLiM"]["generation"] for i in self._tree_sequences]
        #assert len(set(gens)) == 1, ("simulations must be same length (gens).")
        self.generations = gens[0]
        assert self.ancestral_ne, "ancestral_ne not found in metadata."
        assert self.mut, "mut not found in metadata; must enter a mut arg."
        assert self.recomb, "recomb not found in metadata; must enter a recomb arg."

    def _update_tables(self):
        """...Remove extra psuedopopulation nodes."""
        for tree_file in self.trees_files:
            treeseq = pyslim.load(tree_file)

            # get mutable tskit.TableCollection
            tables = treeseq.dump_tables()
            nnodes = tables.nodes.time.size

            # there is a null SLiM population (0) that doesnt really exist
            # and the actual poulation (1). So we set all to 0.
            tables.nodes.population = np.zeros(nnodes, dtype=np.int32)

            # drop nodes that are not connected to anything. This includes
            # the pseudo-nodes representing half of the haploid populations.
            nodes_in_edge_table = list(
                set(tables.edges.parent).union(tables.edges.child)
            )

            # remove the empty population nodes by using simplify, which 
            tables.simplify(
                samples=nodes_in_edge_table,
                keep_input_roots=True,
                filter_individuals=True,
                filter_populations=True,
                filter_sites=False,
            )

            # turn it back into a treesequence
            self._tree_sequences.append(pyslim.load_tables(tables))

    def _match_nodes_and_merge(self):
        """Find matching ancestral nodes between two treeseqs for merging
        ts1 into ts0.

        Given SLiM tree sequences `other` and `ts`, builds a numpy array with length
        `other.num_nodes` in which the indexes represent the node id in `other` and the
        entries represent the equivalent node id in `ts`. If a node in `other` has no
        equivalent in `ts`, then the entry takes the value `tskit.NULL`. The matching
        is done by comparing the IDs assigned by SLiM which are kept in the NodeTable
        metadata. Further, this matching of SLiM IDs is done for times (going 
        backward-in-time) greater than the specified `split_time`.

        Note
        ----
        This must be done before recapitation which would add nodes
        without slim_id metadata.
        """
        ts0 = self._tree_sequences[0]
        ts1 = self._tree_sequences[1]
        split_time = 1000 #int(1 + (2 * self.generations))

        # get an empty array to be filled.
        node_mapping = np.full(ts1.num_nodes, tskit.NULL)

        # get the slim_ids for every node in each ts
        slim_ids0 = np.array([n.metadata["slim_id"] for n in ts0.nodes()])
        slim_ids1 = np.array([n.metadata["slim_id"] for n in ts1.nodes()])

        # get ids of all nodes alive before the split (tskit: think backwards time), 
        # meaning that they are in the ancestor.
        alive_before_split_pop1 = ts1.tables.nodes.time >= split_time

        # get ids of nodes slim_ids for nodes that are in both ts's
        is_1in0 = np.isin(slim_ids1, slim_ids0)

        # get node ids that meet both above booleans
        both = np.logical_and(alive_before_split_pop1, is_1in0)
        sorted_ids0 = np.argsort(slim_ids0)
        matches = np.searchsorted(
            slim_ids0,
            slim_ids1[both],
            side='left',
            sorter=sorted_ids0,
        )
        node_mapping[both] = sorted_ids0[matches]
        match = sum(node_mapping != -1)
        nomatch = sum(node_mapping == -1)
        
        # save it.
        # logger.debug(f"match={match}; nomatch={nomatch}; {self.trees_files}")
        tsu = ts0.union(ts1, node_mapping=node_mapping, check_shared_equality=True)
        self.tree_sequence = pyslim.SlimTreeSequence(tsu)

        # 
        # tables.simplify(
        #     samples=nodes_in_edge_table,
        #     keep_input_roots=True,
        #     filter_individuals=True,
        #     filter_populations=True,
        #     filter_sites=False,
        # )

    def _recapitate(self):
        self.tree_sequence = pyslim.recapitate(
            self.tree_sequence, 
            ancestral_Ne=self.ancestral_ne,
            recombination_rate=self.recomb,
            random_seed=self.seed,
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
            If zero (default) then no subsampling is performed and 
            all nodes are used (in which case replication is not 
            necessary). (Zero subsampling not yet implemented.)
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
        thetas = []

        # get a list of Series
        for rep in range(reps):
            rep_seed = rng.integers(2**31)

            # let toytree do the sampling of N nodes from each pop such that:
            # (pop=pop, time=0, flag=1); this is equivalent to 'alive at' sampling.
            tts = toytree_sequence(self.tree_sequence, sample=sample, seed=rep_seed)
            samples = np.arange(tts.sample[0] + tts.sample[1])
            sample_0 = samples[:tts.sample[0]]
            sample_1 = samples[tts.sample[0]:]
            # logger.debug(f"sample_0 nodes: {sample_0}")
            # logger.debug(f"sample_1 nodes: {sample_1}")

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

            thetas.append(tts.tree_sequence.diversity(sample_0))
            thetas.append(tts.tree_sequence.diversity(sample_1))

        #save all the thetas to a list
        self.thetas = thetas 

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

    def draw_stats(
        self,
        stat: str="divergence",
        window_size: int=20_000,
        sample: Union[int, Iterable[int]]=10,
        reps: int=1,
        seed: Optional[int]=None,
        color: Optional[str]="mediumseagreen",
        ymax: Optional[float]=None,
        plot_style: Optional[str]="fill",
        height: Optional[int]=300,
        width: Optional[int]=500,
        ):
        """Return a toyplot drawing of a statistic across the genome.
        
        If reps > 1 the measurement is repeated on multiple sets of 
        random samples of size `sample`, and the returned statistic
        is the mean with +/- 1 stdev shown. 

        Parameters
        ----------
        """
        # repeat measurement over many random sampled replicates
        rng = np.random.default_rng(seed)
        rep_values = []
        for _ in range(reps):
            rep_seed = rng.integers(2**31)
            tts = toytree_sequence(self.tree_sequence, sample=sample, seed=rep_seed)
            samples = np.arange(tts.sample[0] + tts.sample[1])
            sample_0 = samples[:tts.sample[0]]
            sample_1 = samples[tts.sample[0]:]
            samples = [sample_0, sample_1]

            # select a supported statistic to measure
            if stat == "divergence":
                func = tts.tree_sequence.divergence
            elif stat == "Fst":
                func = tts.tree_sequence.Fst
            else:
                raise NotImplementedError(f"stat {stat} on the TODO list...")
                
            values = func(
                sample_sets=samples,
                windows=np.linspace(
                    start=0, 
                    stop=tts.tree_sequence.sequence_length, 
                    num=round(tts.tree_sequence.sequence_length / window_size)
                )
            )
            rep_values.append(values)

            
        
        # get mean and std
        means = np.array(rep_values).mean(axis=0)
        # stds = np.array(rep_values).mean(axis=0) 
        self.rep_values = rep_values
        self.means = means    

        std_low = []
        std_high = []

        for i in rep_values:
            data = i
            mean_val = np.mean(data)
            low, high = scipy.stats.t.interval(
                alpha=0.95,
                df=len(data) - 1,
                loc=mean_val,
                scale=scipy.stats.sem(data),
            ) 
            std_low.append(low)
            std_high.append(high)

        style = {"fill":str(color)}
        if plot_style == "fill":
            # draw canvas...
            canvas, axes, mark  = toyplot.fill(
                means, height=height, width=width, opacity=0.5, margin=(60, 50, 50, 80),
                ymax=ymax, style=style
            )
        else:
            canvas = toyplot.Canvas(width=width, height=height)
            axes = canvas.cartesian(ymax=ymax, ymin=0)
            for i in range(0, len(means)):
                fill = axes.fill(means[i], means[i]-std_low[i], means[i]-std_high[i], opacity=0.2, style=style)
                line = axes.plot(means, style=style)

        # style axes
        axes.x.ticks.show = True
        axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
        axes.y.ticks.labels.angle = -90
        axes.y.ticks.show = True
        axes.y.ticks.locator = toyplot.locator.Extended(only_inside=True, count=5)        
        axes.label.offset = 20
        axes.label.text = f"{stat} in {int(window_size / 1000)}kb windows"
        return canvas, axes, mark



if __name__ == "__main__":

    import glob
    import shadie
    shadie.set_log_level("DEBUG")

    TREEFILES = sorted(glob.glob(
        "/home/deren/Desktop/standard-params/bryo-mono/"
        "bryo_mono_run1[0-9]_from_smallchrom_2000spo.trees")
    )

    post = OneSim(TREEFILES[0], ancestral_Ne=500, mut=1e-7, recomb=1e-8, chromosome=None)
    print(post.stats())

    post = TwoSims(
        trees_files=[TREEFILES[0], TREEFILES[1]],
        ancestral_ne=200,
        mut=1e-7 / 2.,
        recomb=1e-8,
        chromosome=None,
    )
    print(post.stats(sample=20, reps=20))


    TREEFILES = sorted(glob.glob(
        "/home/deren/Desktop/standard-params/pter-hetero/"
        "pter_hetero_run2[0-9]_from_smallchrom_2000spo.trees")
    )
    post = TwoSims(
        trees_files=[TREEFILES[0], TREEFILES[1]],
        ancestral_ne=200,
        mut=1e-7 / 2.,
        recomb=1e-8,
        chromosome=None,
    )
    print(post.stats(sample=20, reps=20))

    print(post.draw_stats())
