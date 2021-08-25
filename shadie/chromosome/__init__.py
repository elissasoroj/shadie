#!/usr/bin/env python

"""
Utilities for representing a chromosome as a collection of Elements.

The main functions exposed from this subpackage are for generating
Chromosome class objects that can be used to access the Elements as
a DataFrame (:attr:`.data`) and to visualize the genome structure
(:meth:`.inspect` function).

Examples
--------
>>> chrom = shadie.chromosome.default()
>>> print(chrom.data)
>>> chrom.inspect()
"""

from shadie.chromosome.src.factory import random, explicit, default
