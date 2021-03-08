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

You can optionally install SLiM3 [here](https://messerlab.org/slim/). However, for development purposes you can simply refer to the `example_scripts` folder in the `shadie` directory for examples of SLiM scripts, which should be adequate. 


This project has 3 main parts:

1. write a `.slim` script and execute in the program SLiM3 
2. convert a phylogeny (provided as newick string) into SLiM3-comptaible demography 
3. take output from SLiM simulation and provide useful summary statistics, such as MK and dN/dS


#### 1. The .slim Script

See `script.py`

This submodule needs to create a `.slim`  script with correct syntax that will run in SLiM3. Refer to the `README` file in the `/shadie/examples_scripts/` subdirectory for an  example of `.slim script` structure. 


#### 2. Phylogeny --> SliM3 Demography

See `demography.py`

This submodule takes a `newick` file and creates a `.slim` script with population demography that matches the phylogeny. It uses `toytree` to traverse the tree from root to tips and collect information it needs to generate the SLiM demography in a `pandas` dataframe, for example:

| source    | child2    | nodeheight| gen       |
| --------- | --------- | --------- | --------- |
| 0         | 3         | 1000000   | 1         |
| 3         | 5         | 600000    | 801       |
| 0         | 2         | 600000    | 801       |
| 3         | 4         | 300000    | 1401      |
| 0         | 1         | 300000    | 1401      |
 
- each node moving from tips to root is renamed with the lowest number of the child tips
- `gen` = 1+(abs(nodeheight-treeheight)*(gentime/treeheight))
- This example data is from a tree generated in toytree, which has no outgroup. Because the first split happens at gen1, we probably want a burn-in time before the simulation begins. 

 This will be used to generate this part of the script, which controls when populations in the program split into subpopulations:

```
 #write beginning row:
1 { sim.addSubpop("p1", K); 
}
#for row in pandas DF, write the following:
{gen}{
    sim.addSubpopSplit("p{dest}", K, p{source});
    p{source}.setMigrationRates(p{dest}, 1.0);
}
{gen+1}{
    p1.setMigrationRates(p2, 0.0);
}
#append until out of rows

#then, write last line:
self.gentime late() { 
    sim.outputFull(); 
}
```


#### 3. Summary Statistics

*Not yet written*



