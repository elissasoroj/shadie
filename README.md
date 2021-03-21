# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around SLiM3 that implements selection on alternating hapliod/diploid lifecycles
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography

### In Development

To install locally:
```bash
#dependencies:
conda install [pandas, numpy, toyplot, toytree, altair, msprime, pyslim, loguru] -c conda-forge

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

Right now, we are focusing on **#1 and #2**. See detail below ~

---
### 1. Phylogeny --> SliM3 Demography

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
~

**Working Example**

Right now, the code can successfully traverse a tree and save the output to a pandas DF:

1. install `shadie` from the `./shadie` local directory
```bash
pip install -e .
```

2. launch jupyter notebook and create a new Python3 Notebook
```bash
jupyter notebook
```

3. In the notebook:
```python
#import dependencies
from shadie import Demography
import toytree
import numpy as np

#create a test tree using toytree
tree = toytree.rtree.unittree(ntips=10, treeheight=1e6, seed=123)

# set Ne values on the nodes of the tree
tree = tree.set_node_values(
    feature="Ne", 
    values={i: np.random.randint(10000, 100000) for i in tree.idx_dict},
)

#retrieve demography from the tree and save as pandas DF:
dem = Demography(tree, Ne=10000)
Demography.get_demog(dem)
```
~

**To Add**
1. Read 'Ne' from each tree node and store in `self.demog`
2. Store all the SLiM lines in a dictionary 
	a. `(key, value) = (generation, (str))` where (str) = lines of code that need to be printed to the SLiM script
---

### 2. Chromosome

See `chromosome.py` for script

**Working Example:**

1. install `shadie` from the `./shadie` local directory
```bash
pip install -e .
```

2. launch jupyter notebook and create a new Python3 Notebook
```bash
jupyter notebook
```

3. In the notebook:
```python
from shadie import Chromosome

#create a Chromosome class object
#default genome size = 1e6; set size with "genome_size=(int)"
random_chromosome = Chromosome()

#build a random chromosome with `make_rand` function
random_chromosome.make_rand()

#inspect your chromosome using 'review' function
#dataframe of genomic elements
random_chromosome.review("elements")

#Chromosome plot with # of genes, mean intron & exon length, etc...
random_chromosome.review("chromosome")

#interactive `altair` plot:
random_chromosome.review("interactive")
```
~

**To Add**
1. Interactive `altair` chart
	
	a. Tooltips: hovering over chart displays:
	* ~~genomic element type (e.g. coding, non-coding, etc...)~~
	* ~~start and stop bases~~
	* element length
	* gene #
	* mutation types 

	b. double-plot zxoom feature
	* one full chromosome plot with selector rectangle (see the last section of [this chapter of altair docs](https://altair-viz.github.io/user_guide/interactions.html#conditions-making-the-chart-respond))
	* selector defines coordinate for the plot below, which displays a zoomed in version of the chromosome

	c. choose better colors

2. Add more customization to `make_rand()`

	a. # of genes

	b. average # of exons and introns per chromosome (adjusts chance generator adds another non-coding regiono)

	c. average length of exons and introns (adjusts distributions that exon and intron lengths are drawn from)

	d. average length of non-coding regions (adjusts distirbution that non-coding regions are drawn from)

	e. Need to validate default values and provide average numbers (with citation) for users to reference

3. Load chromosome structure from user dictionary 

	a. user can provide a chromosome structure as dictionaries types directly into the API

	- dictionary converted to pandas df

	b. they can also provide a dataframe or .csv file

	c. `Chromosome.make()` reads the dataframe and generate:
	- mutation types
	- genomic element types
	- genomic elements (chromosome structure)

---
### 3. Generate the .slim Script

See `shadie.py` for script

This submodule needs to create a `.slim`  script with correct syntax that will run in SLiM3. Refer to the `README` file in the `/shadie/examples_scripts/` subdirectory for an  example of `.slim script` structure. 

**Steps**
1. Read in Chromosome object
	
	a. mutation rate

	b. mutation types

	c. genomic element types

	d. genomic elements (chromosome)

2. Read in demography dictionary

3. User will choose `organism` type, which determines reproduction

4. `Shadie()` uses information from `Chromosome` to write initialize callbacks to a .slim file

5. `Shadie()` uses `Demography` dictionary, pre-defined `reproduction` dictionaries, and `early/late` callback dictionaries to construct the rest of the .slim script

	a. possible future functionality to add: let users adjust reproduction settings

6. Shadie() passes the script to SliM3 and runs the simulation

7. Shadie() reads the output back in, allowing the user to look at output with `Postsim` module

---
### 4. Summary Statistics (Postsim module)

*Not yet written*
