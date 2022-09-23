## First...

a. make sure you have SliM3, Python, and Jupyter notebooks installed on your machine.

b. Deploy a Jupyter in terminal and create a new notebook.

```bash
jupyter notebook
```

## 1. Make Custom Mutations

The first step in setting up a SLiM simulation is to define the mutations that can occur on the chromosome of the organisms in your simulation. Each mutation will have a **dominance coefficient** (h) and a **fitness effect** (s). The dominance coefficient will determine how the fitness of a heterozygote is calculated. 

The fitness effect is a `float` value between 0 and 1 inclusive. In SLiM, the multiplicative fitness effect of a mutation with selection coefficient `s` is `1+s` for a homozygote; in a heterozygote it is `1+hs`, where `h` is the dominance coefficient. 

`shadie` automatically sets all neutral mutations to be 'ignored' by the SLiM3 simulation. `shadie` will use `pyslim` after the SLiM simulation is complete to overlay neutral mutations using a much faster coalescent method. 

## 2. Create Genomic Elements

## 3. Build Chromosome

## 4. 