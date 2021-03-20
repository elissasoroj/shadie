# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around SLiM3 that implements selection on alternating hapliod/diploid lifecycles
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography

### In Development

To install locally:
```bash
#dependencies:
conda install [pandas, numpy, toyplot, toytree, msprime, pyslim, loguru] -c conda-forge

#clone and install
git clone [https://github.com/elissasoroj/shadie.git]
cd ./shadie
pip install -e .
```

You can optionally install SLiM3 [here](https://messerlab.org/slim/). However, for development purposes you can simply refer to the `example_scripts` folder in the `shadie` directory for examples of SLiM scripts, which should be adequate. *You do NOT have to install SliM3*


This project has 3 main parts:

1. convert a phylogeny (provided as newick string) into SLiM3-comptaible demography 
2. generate chromosome, either random or based on user input
3. write a `.slim` script and execute in the program SLiM3 
4. take output from SLiM simulation and provide useful summary statistics, such as MK and dN/dS

Right now, we are focusing on **#2**. See detail below ~


#### 1. Phylogeny --> SliM3 Demography

See `demography.py` for script

This submodule takes a `Toytree` tree object and creates a `.slim` script with population demography that matches the phylogeny. It uses `Toytree` to traverse the tree from root to tips and collect information it needs to generate the SLiM demography in a `pandas` dataframe, for example:

 |   gen   | source  | child0   | child1   | Ne      |
 |-------- | ------- | -------- | -------- | ------- |
 |1        | 9       | 8        | 7        | 1000    |
 |501      | 7       | 6        | 5        | 2000    |
 |901      | 8       | 4        | 3        | 2000    |
 |1401     | 5       | 2        | 1        | 4000    |
 
- each node is named by Toytree from tips to root; `Demography` retains that naming
- `gen` = int((theight - node.height) + 1)

 This will be used to generate this part of the script, which controls when populations in the program split into subpopulations:

```java
//write beginning row:
1 { sim.addSubpop("{root}", {Neroot}); }

// for row in pandas DF, write the following:
{gen} { 
	sim.addSubpopSplit("{child0}", {Nechild0}, p1);
    sim.addSubpopSplit("{child1}", {Nechild1}, p1);}
{gen}:{gen+1} {
    {root}.setMigrationRates(c({child0}, {child1}), c(0.5, 0.5));} //may not be necessary
{gen+1} {
    {root}setSubpopulationSize(0);}
// append until out of rows

//then, write last line:
10000 late() { sim.outputFull(); }
```


#### 2. Chromosome

See `chromosome.py` for script

Working Example:

install `shadie` from the `./shadie` local directory
```bash
pip install -e .
```

launch jupyter notebook and create a new Python3 Notebook
```bash
jupyter notebook
```

In the notebook:
```python
from shadie import Chromosome

#create a Chromosome class object
#default genome size = 1e6; set size with "genome_size=(int)"
random_chromosome = Chromosome()

#build a random chromosome with `make_rand` function
Chromosome.make_rand(random_chromosome)

#inspect your chromosome using 'review' function
#dataframe of genomic elements
Chromosome.review(random_chromosome, item = "elements")

#Chromosome plot with # of genes, mean intron & exon length, etc...
Chromosome.review(random_chromosome, item = "chromosome")
```


#### 3. Generate the .slim Script

See `shadie.py` for script

This submodule needs to create a `.slim`  script with correct syntax that will run in SLiM3. Refer to the `README` file in the `/shadie/examples_scripts/` subdirectory for an  example of `.slim script` structure. 



#### 4. Summary Statistics

*Not yet written*
