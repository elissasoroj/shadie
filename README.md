# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around **SLiM3** that provides a user-friendly python API interface, implements selection on alternating hapliod/diploid lifecycles.

A newer version of SLiM, SLiM4 has been released. `shadie` is being updated for compatibility with SLiM4. For now, if you want to use `shadie` you should wait to upgrade. 

Docs (under construction): https://elissasoroj.github.io/shadie/


### *`shadie` is Under Active Development*

To install locally:
```bash
#dependencies:
conda install [pandas, numpy, toyplot, toytree, altair, msprime, pyslim, tskit, loguru, IPython] -c conda-forge

#clone and install from GitHub
git clone [https://github.com/elissasoroj/shadie.git]
cd ./shadie
pip install -e .
```

You can install SLiM3 [here](https://messerlab.org/slim/). *You must install SliM3 for most of the functionality of `shadie` to work.* Note: `shadie` does *not* yet work with the newest release, SLiM4. 

~

**Dependencies:**

* numpy
* pandas
* toytree
* toyplot
* altair
* pyslim
* tskit

~

Planned implementation:
* convert user-provided phylogeny into SLiM3-compatible subpopulation demography. 
* post-sim analysis module