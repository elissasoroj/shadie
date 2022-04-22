# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around SLiM3 that provides a user-friendly python API interface, implements selection on alternating hapliod/diploid lifecycles,
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography. 

Documentation is forthcoming.

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

You can install SLiM3 [here](https://messerlab.org/slim/). For development purposes you can refer to the `example_scripts` folder in the `shadie` directory for examples of SLiM scripts. It is  highly recommended you also familiarize yourself with the SliM language `Eidos`; the SLiMgui is included with the Mac OS install and may be of tremendous help in this regard.  *You must install SliM3 for most of the functionality of `shadie` to work*

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
