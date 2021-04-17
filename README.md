# SHaDiE (Simulating Haploid-Diploid Evolution)

`shadie` is a wrapper around SLiM3 that provides a user-friendly python API interface, implements selection on alternating hapliod/diploid lifecycles,
and converts user-provided phylogeny into SLiM3-compatible subpopulation demography

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

**This project has 5 main parts:**

1. convert a phylogeny (provided as newick string) into SLiM3-comptaible demography 
2. generate chromosome, either random or based on user input
3. pre-define reproduction scripts that reflect lifecycles of major plant lineages
4. write a `.slim` script and execute in the program SLiM3 
5. take output from SLiM simulation and provide useful summary statistics, such as MK and dN/dS

Right now, we are focusing on **#4 and #5**. See detail below ~

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

 This will be used to generate this part of the script, which controls when populations in the program split into subpopulations.

~

**Working Example**

*Refer to notebook #3 in the `notebooks` directory.*

Right now, the code can successfully traverse a tree and save the output to a pandas DF:

1. launch jupyter notebook and create a new Python3 Notebook
```bash
jupyter notebook
```

2. In the notebook:
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
dem = Demography(tree)
Demography.get_demog(dem)
```

---

### 2. Chromosome

See `chromosome.py` for script. 

*Refer to notebook #2 in the `notebooks` directory.*

**Working Example:**

In a Jupyter Notebook:

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
	* ~~element length~~
	* gene #
	* mutation types 

	~~b. double-plot zoom feature~~
	* one full chromosome plot with selector rectangle (see the last section of [this chapter of altair docs](https://altair-viz.github.io/user_guide/interactions.html#conditions-making-the-chart-respond))
	* selector defines coordinate for the plot below, which displays a zoomed in version of the chromosome

	~~c. choose better colors~~

2. Add more customization to `make_rand()`

	~~a. # of genes~~

	b. average # of exons and introns per chromosome (adjusts chance generator adds another non-coding region)

	c. average length of exons and introns (adjusts distributions that exon and intron lengths are drawn from)

	d. average length of non-coding regions (adjusts distirbution that non-coding regions are drawn from)

	~~e. Need to validate default values and provide average numbers (with citation) for users to reference~~

3. Load chromosome structure from user dictionary 

	a. user can provide a chromosome structure as dictionaries types directly into the API

	- dictionary converted to pandas df

	b. they can also provide a dataframe or .csv file

	c. `Chromosome.make()` reads the dataframe and generate:
	- mutation types
	- genomic element types
	- genomic elements (chromosome structure)

---
### 3. Reproduction 

*Not yet written*

See `reproduction.py` for script

This submodule will pre-define script parameters for different life-cycles. 

For now, reproduction is an optional argument in the `Shadie()` class object that will add alternation of generations to the simulation:

~

**Working Example**

```python
#imports
from shadie import Shadie

#init a Shadie class object with reproduction settings
altgen_sim = Shadie(reproduction = True, model = "nonWF")

altgen_sim.write()
altgen_sim.run()

#see the next section for more details on write/run

```

---
### 4. Generate the .slim Script

See `main.py` for script

This submodule needs to create a `.slim`  script with correct syntax that will run in SLiM3. Refer to the `README` file in the `/shadie/examples_scripts/` subdirectory for an  example of `.slim script` structure. 

**Steps**
~~1. Read in Chromosome object~~
	
	a. mutation rate

	b. mutation types

	c. genomic element types

	d. genomic elements (chromosome)

2. Read in demography dictionary

	This part of the script works, but is not compatible with the `nonWF` model required by complex `reproduction` required to model plant  lifecycles. We are working on a different workflow, tentatively outlines below:

	a. `shadie` uses the demography df to set up parallel simulations, one for each edge of the tree. 

	b. Ideally, these simulations will be able to run in  parallel 

	c. *Note* this will be a separate program

3. User will choose `organism` type, which determines reproduction

4. ~~`Shadie()` uses information from `Chromosome` to write initialize callbacks to a .slim file~~

5. `Shadie()` uses `Demography` dictionary, pre-defined `reproduction` dictionaries, and `early/late` callback dictionaries to construct the rest of the .slim script

	a. possible future functionality to add: let users adjust reproduction settings rather than just choosing from defaults

6. ~~Shadie() passes the script to SliM3 and runs the simulation~~

7. ~~Shadie() reads the output back in, allowing the user to look at output with `Postsim` module~~

**Future Changes**:
* Remove nucleotides from the simulation (should run faster without them and they are no longer needed)
* Re-configure reproduction callbacks to write script using `Reproduction` class instead of `Shadie`
* Make all scripts consistent with nonWF models
	* consider making a separate module for writing simple WF models?
	* *or* just implement WF into nonWF model

~

**Working Example:**
*Refer to notebook #4 in the `notebooks` directory.*

```python
from shadie import Chromosome

sim1 = Shadie() #this initializes a simulation with all the shadie defaults

#write a shadie.slim file to the working dir
sim1.write()

#run the shadie.slim file
sim1.run
```

---

### 5. Summary Statistics (`Postsim` module)

This module will read the output of SLiM (reads in tree sequence) and overlays neutral mutations uing `pyslim`. 
It then allows the user to analyze their simulation results. 

**To Add:**

1. View phylogeny with tip `idx` that allows users to refer to different lineages

2. `sample()` function: will sample individuals from the tips (last generation of the simulation)

	a. random sample

	b. lineage-specific sample

	c. more control over sample (e.g. plot a distribution of individuals based on genetic variation from ancestral state and choose based on degree of variance)

3. dNdS calculation

4. Mk-test (needs phylogeny)

5. Linkage Disequilibrium

6. ~~Plot mutations on chromosome map to show postitions/densities~~

7. Allow users to read in a `.trees` file plus saved info from the simulation
	* alternatively, we couls have `shadie` "read" the .slim file and save all the info to new objects instead

~

**Working Example:**
*Refer to notebook #5 in the `notebooks` directory.*

```python
from shadie import PostSim

sim1_post = PostSim(sim1)

#see a plot of all the mutations that occurred (this plots ALL the mutations)
sim1_post.summary()

#calculate non-synonymous:synonymous mutation ratio (will become dNdS later)
sim1_post.dNdS()
```
