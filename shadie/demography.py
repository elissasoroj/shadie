#!/usr/bin/env python

"""
Traverses phylogeny from tips to root and builds corresponding 
population demography for SLiM
"""

#imports
import pandas as pd
import numpy as np
import toytree
from loguru import logger


class Demography:
    """
   
    Parameters
    ----------
    tree: (str, Toytree)
        Ultrametric ToyTree tree object with units in generations.
    """
    def __init__(self, tree, Ne=None):
       
        # store input params
        self.tree = toytree.tree(tree)

        # apply Ne
        self.tree = self.tree.set_node_values("Ne", default=Ne)
       
        # will be used to store output results
        self.demog = None #demography in proper format


    def get_demog(self):
        """
        Returns demography scenario based on an input tree and admixture
        edge list with events in the format (source, dest, start, end,
        rate). Time on the tree is defined in units of generations.

        demogdf = \
           gen    src  child0   child1   Ne     
        0   100    1      2         3     xxx
        1   1000   2      4         5     xxx
        2   1200   3      6         7     xxx

        """
        # a dataframe for storing info above.
        demogdf = pd.DataFrame(
            columns=['gen', 'src', 'child0', 'child1', 'Ne'],
            data=None,
        )

        # total tree height
        theight = self.tree.treenode.height

        # traverse tree from root to tips
        for idx, node in enumerate(self.tree.treenode.traverse('preorder')):
            if len(node.children) == 2: #patrick told me this avoids silly errors
                gen = int((theight - node.height) + 1)
                demogdf.loc[idx, 'gen'] = gen
                demogdf.loc[idx, 'src'] = f'p{node.idx}'
                demogdf.loc[idx, 'child0'] = f'p{node.children[0].idx}'
                demogdf.loc[idx, 'child1'] = f'p{node.children[1].idx}'
                demogdf.loc[idx, 'Ne'] = node.Ne
                #I started by following the f-string format, but it seems like just having the value is better
        self.demog = demogdf
        logger.debug(self.demog)


#This is the .slim script structure that needs to be written by the submodule:
   #Plan is to reference pandas DF with columns: 
    #gen, src, child0, child1, Ne
    # sorted by youngest --> oldest gen (or longest --> shortest nodeheight)

"""
#write beginning row:
1 { sim.addSubpop("{root}", {Neroot}); }

#for row in pandas DF, write the following:
{gen} { sim.addSubpopSplit("{child0}", {Nechild0}, p1);
    sim.addSubpopSplit("{child1}", {Nechild1}, p1);}
{gen}:{gen+1} {
    {root}.setMigrationRates(c({child0}, {child1}), c(0.5, 0.5));
}
{gen+1} {
    {root}setSubpopulationSize(0);
}
#append until out of rows

#then, write last line:
10000 late() { sim.outputFull(); }
"""


if __name__ == "__main__":

    # make a random tree with 10 tips and root height 1M
    tree = toytree.rtree.unittree(ntips=10, treeheight=1e6, seed=123)

    # set Ne values on the nodes of hte tree
    tree2 = tree.set_node_values(
        feature="Ne", 
        values={i: np.random.randint(10000, 100000) for i in tree.idx_dict},
    )

   # TODO: dem loads a tree and parses Ne from nodes
    dem2 = Demography(tree2)
    Demography.get_demog(dem2)

    # supported already
    dem = Demography(tree, Ne=10000)
    Demography.get_demog(dem)
  