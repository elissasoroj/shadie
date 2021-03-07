#!/usr/bin/env python

"""
uses msprime to run neutral burn-in then saves to .trees to be used as starting state for SLiM simulation
imports .trees from SLiM simulation into msprime to overlay neutral mutations using coalescent

"""

#imports
import msprime
import pyslim
