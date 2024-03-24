# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a Python wrapper around the evolutionary simulation program 
**SLiM** that provides a user-friendly API interface. In addition to 
implementing standard SLiM operators, `shadie` includes a rich library of
functions to generate genome structure (chromosomes), mutation types, and
selection scenarios. In particular, this includes the `shadie.reproduction`
module for creating simulations of alternating hapliod/diploid lifecycles.

`shadie` has been updated for compatibility with SLiM4.

Docs (under construction): https://elissasoroj.github.io/shadie/

### *`shadie` is Under Active Development*

As such, the current release of `shadie` may not always be working. Please 
don't hesitate to start a discussion if you have any questions. 


### Installation

To install locally:
```bash
# recommended installation method (installs shadie and slim>=4.2)
conda install shadie -c conda-forge
```

Alternatively, for developers:
```bash
conda install shadie -c conda-forge --only-deps
git clone https://github.com/elissasoroj/shadie.git
cd ./shadie/
pip install -e . --no-deps
```

<!-- You can install SLiM4 [here](https://messerlab.org/slim/). `shadie` now requires the latest release of SLiM 4.1.  -->

~

**Dependencies:**

* numpy>=1.26.2
* pandas>=2.1.4
* toytree>=3.0.0
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
