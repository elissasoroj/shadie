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


class PureSlim:
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
        ancestral_burnin: str,
        seed: Optional[int]=None,
        **kwargs,
        ):

        # hidden attributes
        self._tree_files: List[str] = tree_files
        self._tree_sequences = [pyslim.load(i) for i in self._tree_files]
        self._nts: int = len(self._tree_files)

        self._burnin= pyslim.load(ancestral_burnin)

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
        self._union_all()

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

    def _union_all(self):
        """Calls tskit union, simplify, mutate to complete neutral sim.
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

        self.tree_sequence = tsu

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
	    