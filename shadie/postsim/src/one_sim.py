#!/usr/bin/env python

"""
A returned object class from a shadie simulation call.
"""

from typing import List, Optional, Union, Iterable
import numpy as np
import pyslim
import tskit
import msprime
import toyplot
import pandas as pd
import numpy as np
import scipy.stats
from loguru import logger

from toytree.utils.tree_sequence import ToyTreeSequence
from toytree.utils.utils import ScrollableCanvas
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
        tree_sequence: 'tskit.trees.TreeSequence',
        seed: Optional[int]=None,
        **kwargs,
        ):

        # hidden attributes
        self.tree_sequence = pyslim.load(tree_sequence)
        """A SlimTreeSequence that has been recapitated and mutated."""

        # attributes to be parsed from the slim metadata
        self.generations: int=0
        """The SLiM simulated length of time in diploid generations."""
        self.ancestral_Ne: int=kwargs.get("ancestral_Ne")
        """The SLiM simulated diploid carrying capacity"""
        self.recomb: float=kwargs.get("recomb")
        """The recombination rate to use for recapitation."""
        self.mut: float=kwargs.get("mut")
        """The mutation rate to use for recapitated treesequence."""
        self.chromosome: ChromosomeBase=kwargs.get("chromosome")
        """The shadie.Chromosome class representing the SLiM genome."""
        self.rng: np.random.Generator=np.random.default_rng(seed)

        # try to fill attributes by extracting metadata from the tree files.
        self._extract_metadata()
        self._tskit_complete()

    def _extract_metadata(self):
        """Extract self attributes from shadie .trees file metadata.

        TODO: can more of this be saved in SLiM metadata?
        """
        self.generations = self.tree_sequence.metadata["SLiM"]["generation"]
        assert self.ancestral_Ne, "ancestral_Ne not found in metadata; must enter an ancestral_Ne arg."
        assert self.mut, "mut not found in metadata; must enter a mut arg."
        assert self.recomb, "recomb not found in metadata; must enter a recomb arg."

    def _tskit_complete(self):
        """Calls tskit union, simplify, mutate to complete neutral sim.
        """
        self._update_tables()
        self._recapitate()
        self._mutate()

    def _update_tables(self):
        """DEPRECATED.

        Alternative approach to remove nulll pop and divide time.
        This didn't work, still squishes edges...

        Try this stuff next...
        https://github.com/tskit-dev/pyslim/blob/625295ba6b4ae8e8400953be65b03b3630c1430f/docs/vignette_continuing.md#continuing-the-simulation
        """
        # get mutable tskit.TableCollection
        tables = self.tree_sequence.dump_tables()
        nnodes = tables.nodes.time.size

        # there is a null SLiM population (0) that doesnt really exist
        # and the actual poulation (1). So we set all to 0.
        tables.nodes.population = np.zeros(nnodes, dtype=np.int32)
        meta = tables.metadata
        meta["SLiM"]["generation"] = int(meta["SLiM"]["generation"] / 2.)
        tables.metadata = meta
        tables.nodes.time /= 2
        tables.mutations.time /= 2.

        tables.simplify(
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
            tts = ToyTreeSequence(self.tree_sequence, sample=sample, seed=seed)
            samples = np.arange(tts.sample[0] + tts.sample[1])
            sample_0 = samples[:tts.sample[0]]
            sample_1 = samples[tts.sample[0]:]

            stats = pd.Series(
                index=["theta_0", "theta_1", "Fst_01", "Dist_01", "D_Taj_0", "D_Taj_1"],
                name=str(rep),
                data=[
                    tts.tree_sequence.diversity(sample_0_),
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
