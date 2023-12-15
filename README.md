# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around **SLiM3** that provides a user-friendly python API interface, implements selection on alternating hapliod/diploid lifecycles.

Docs (under construction): https://elissasoroj.github.io/shadie/


### *`shadie` is Under Active Development*

As such, the current release of `shadie` may not always be working. Please don't hesitate to start a discussion if you have any questions. 


### Installation

To install locally:
```bash
#dependencies:
conda install [pandas, numpy, toyplot, toytree, altair, msprime, pyslim, tskit, loguru, scipy] -c conda-forge

#clone and install from GitHub
git clone https://github.com/elissasoroj/shadie.git
cd ./shadie
pip install -e .

#shadie requires utilities from the development version of the package `toytree`
#clone and install from GitHub
git clone https://github.com/eaton-lab/toytree.git
cd ./toytree
git checkout toy3
pip install -e .
```

You can install SLiM4 [here](https://messerlab.org/slim/). `shadie` now requires the latest release of SLiM 4.1. 

~

**Dependencies:**

* numpy>=1.26.2
* pandas>=2.1.4
* toytree - dev branch toy3
* toyplot>=1.0.3
* altair>=5.2.0
* pyslim>=1.0.2
* msprime>=1.3.0
* tskit>=0.5.6
* loguru>=0.6.0

~

Planned implementation:
* additional models (currently implementing: red algae, *Vittaria appalachiana*)
* convert user-provided phylogeny into SLiM3-compatible subpopulation demography. 
* post-sim analysis module
