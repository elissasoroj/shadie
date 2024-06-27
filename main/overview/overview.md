<h1>Welcome to the SHaDiE Documentation</h1>

## Overview

`shadie` (Simulating Haploid-Diploid Evolution) is a Python API wrapper for the forward-in-time evolutionary simulation software, SLiM. `shadie` simulates non-neutral mutations in SLiM, then uses functionality provided by `pyslim` and `tskit` to overlay neutral mutations using coalescent simulations in `msprime`. 

`shadie` also provides a user-friendly Python API framework for using SLiM. This is meant to make SLiM a little more approachable for users with prior Python knowledge, however it does not incorporate all the functionality of SLiM and it is hightly encouraged for users to familiarize themselves with the SLiMgui and the SliM Manual. The [Messer Lab](https://messerlab.org/slim/) provides excellent documentation for learning SLiM. 

Finally, `shadie` was first developed with the intention of modeling realistic plant lifecycles in SLiM, so much of the tutorials and documentation focus on this use case. 

---

## Use

`shadie` works best in [Jupyter Notebooks](https://jupyter.org/). There are plans to develop a command line tool, however `shadie` relies heavily on visualizations to aid users in designing their SLiM model and analyzing the output of their simulations. For this reason a Python API is recommended, especially for new users. 

### SLiM

[SLiM](https://messerlab.org/slim/) is powerful evolutionary simulation framework developed by the Messer Lab. They provide fanastic documentation and resources on their website. 

### pyslim

[`pyslim`](https://pyslim.readthedocs.io/en/latest/index.html) is a package that modifies `tskit` tree sequence files outpu by SLiM. `shadie` utilizes `pyslim` to overlay neutral muations with `msprime`. This method uses coalescent and is therefore *much* faster than SLiM in most cases. 

### Dependencies

**General**

	- subprocess
	- numpy
	- pandas
	- SLiM
	- msprime
	- pyslim
	- tskit
	- loguru

**Plotting:**

	- toyplot: for static plots
	- toytree: for static plots
	- altair: for interactive plots
	- notebook: for displaying plots in Jupyter Notebooks
