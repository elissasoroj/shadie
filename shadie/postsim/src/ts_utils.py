#!/usr/bin/env python

"""
A returned object class from a shadie simulation call.
"""

from typing import Optional, Union, Iterable, List
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
import pyslim
import tskit
import msprime
import toyplot
import tskit
import scipy.stats
from loguru import logger

from toytree.utils.src.toytree_sequence import ToyTreeSequence
from shadie.chromosome.src.classes import ChromosomeBase

logger = logger.bind(name='shadie')


def stats(
        tree_sequence,
        sample: int = 10,
        seed: Optional[int]=None,
        reps: int=100
        ):
        """Calculate statistics summary on pure TreeSequence.
        
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
            tts = toytree_sequence(tree_sequence, sample=sample, seed=seed)
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

def draw_stats(
        tree_sequence,
        stat: str="diversity",
        window_size: int=500,
        sample=15,
        reps: int=200,
        seed=None,
        color="lightseagreen"
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
            func = tree_sequence.diversity
        else:
            raise NotImplementedError(f"stat {stat} on the TODO list...")

        # repeat measurement over many random sampled replicates
        rng = np.random.default_rng(seed)
        rep_values = []
        for _ in range(reps):
            ndt = tree_sequence.tables.nodes
            mask = (ndt.population == 0) & (ndt.time == 0) & (ndt.flags == 1)
            arr = np.arange(mask.shape[0])[mask]
            size = min(arr.size, sample)
            samples = rng.choice(arr, size=size, replace=False)

            values = func(
                sample_sets=samples,
                windows=np.linspace(
                    start=0, 
                    stop=tree_sequence.sequence_length, 
                    num=round(tree_sequence.sequence_length / window_size)
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
        axes.y.ticks.locator = toyplot.locator.Extended(only_inside=True, count=8)        
        axes.label.offset = 20
        axes.label.text = f"{stat} in {int(window_size)}bp windows"
        return canvas, axes, mark