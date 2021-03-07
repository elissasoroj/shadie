# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around SLiM3 that implements selection on alternating hapliod/diploid lifecycles
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography

### In Development

To install locally:
```bash
#dependencies:
conda install [pandas, numpy, toytree, msprime, pyslim] -c conda-forge

#clone and install
git clone [https://github.com/elissasoroj/shadie.git]
cd ./shadie
pip install -e .
```

You can optionally install SLiM3 (here)[]. However, for development purposes you can simply refer to the `example_scripts` folder in the `shadie` directory for examples of SLiM scripts, which should be adequate. 


This project has 3 main parts:

1. write a `.slim` script and execute in the program SLiM3 
2. convert a phylogeny (provided as newick string) into SLiM3-comptaible demography 
3. take output from SLiM simulation and provide useful summary statistics, such as MK and dN/dS


#### 1. The .slim Script

This submodule needs to create a `.slim`  script with correct syntax that will run in SLiM3. Refer to the `README` file in the `/shadie/examples_scripts/` subdirectory for an  example of `.slim script` structure. 

The 