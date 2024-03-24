#!/usr/bin/env python

"""Code to generate Reproduction blocks in shadie code given functions
and parameter arguments specific to organismal clades.


SLiM population codes
======================
Codes used to refer to individuals from specific populations and/or
life stages in different organismal simulations. These encodings are
defined in: bryo_scripts.py, fern_scripts.py, etc.

Bryo dioicous Bryophytes
-----------------------
p0 = haploid population
p1 = diploid population
tag=1 = gametophyte (1N)
tag=2 = gametophyte clones (1N)
tag=3 = sporophyte (2N)
tag=4 = sporophyte clones (2N) -- None in this model
tagL0 = Female = True, Male = False

Bryo monoicous Bryophytes
-----------------------
p0 = haploid population
p1 = diploid population
tag=1 = gametophyte (1N)
tag=2 = gametophyte clones (1N)
tag=3 = sporophyte (2N)
tag=4 = sporophyte clones (2N) -- None in this model
tagL0 = Female = True, Male = False

Fern pteridophyte homosporous
-----------------------------
p0 = haploid population
p1 = diploid population
tag=1 = gametophyte (1N)
tag=2 = gametophyte clones (1N)
tag=3 = sporophyte (2N)
tag=4 = sporophyte clones (2N)
tagL0 = Female = True, Male = False

Fern pteridophyte heterosporous
-----------------------------
p0 = haploid population
p1 = diploid population
tag=1 = gametophyte (1N)
tag=2 = gametophyte clones (1N)
tag=3 = sporophyte (2N)
tag=4 = sporophyte clones (2N)
tagL0 = Female = True, Male = False

Angio angiosperm
-----------------------------
p0 = haploid population
p1 = diploid population
tag=>1M- tmp parental tags
tag=<2,000,000 = female gametophyte (1N)
tag=>2,000,000 = male gametophyte (1N)
tag2 = gametophyte clones (1N) tag -------------- None in this model
tag3 = sporophyte (2N) tag
tag4 = sporophyte clones (2N) tag
tag0 = used pollen (1N) tag
"""
