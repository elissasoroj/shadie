<h1>Installation</h1>

## Local Install

The recommended way to install shadie is to use `conda`. This will ensure that
all dependencies, including an up-to-date version of SLiM, are also installed.

```bash
# recommended installation method (installs shadie and slim>=4.2)
conda install shadie -c conda-forge
```

Alternatively, for developers:
```bash
conda install shadie -c conda-forge --only-deps
git clone https://github.com/elissasoroj/shadie.git
cd ./shadie/
pip install -e . --no-deps
```
<!-- #dependencies:
conda install pandas, numpy, toyplot, toytree, altair, msprime, pyslim, tskit, \
	loguru, ipython -c conda-forge

#clone and install from GitHub
git clone [https://github.com/elissasoroj/shadie.git]
cd ./shadie
pip install -e .
``` 
-->

## Troubleshooting

If you need help, please use the [Discussions](https://github.com/elissasoroj/shadie/discussions) on GitHub. If you want to request a feature or report a bug, please open an [issue](https://github.com/elissasoroj/shadie/issues).
