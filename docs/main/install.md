<h1>Installation</h1>

## Local Install

Currently, the only way to install `shadie` is through a local install:

```
#dependencies:
conda install [pandas, numpy, toyplot, toytree, altair, msprime, pyslim, tskit, loguru, scipy] -c conda-forge

#clone and install from GitHub
git clone https://github.com/elissasoroj/shadie.git
cd ./shadie
pip install -e .

#shadie requires utilities from the development version of the package `toytree`
#clone and install from GitHub
git clone https://github.com/eaton-lab/toytree.git
cd ./toytree
git checkout toy3
pip install -e .
```

## Troubleshooting

If you need help, please use the [Discussions](https://github.com/elissasoroj/shadie/discussions) on GitHub. If you want to request a feature or report a bug, please open an [issue](https://github.com/elissasoroj/shadie/issues).
