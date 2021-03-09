#!/usr/bin/env python

"""
Traverses phylogeny from tips to root and builds corresponding 
population demography for SLiM
"""

#imports
import toytree
import pandas as pd


class Demography:
    """
    
    Parameters
    ----------
    tree: (str, Toytree)
        Ultrametric ToyTree tree object with units in generations.
    """
    def __init__(self, tree, Ne=1000):
        
        # store input params
        self.tree = toytree.tree(tree)

        # apply Ne
        self.tree = self.tree.set_node_values("Ne", default=Ne)
        
        # will be used to store output results
        self.demog = None #demography in proper format


    def _get_demog(self):
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
            if node.children:
                gen = int((theight - node.height) + 1)
                demogdf.loc[idx, 'gen'] = gen
                demogdf.loc[idx, 'src'] = f'p{node.idx}'
                demogdf.loc[idx, 'child0'] = f'p{node.children[0].idx}'
                demogdf.loc[idx, 'child1'] = f'p{node.children[1].idx}'
        self.demog = demogdf


    def demog_to_slim(self):
        """
        
        """

        
        # # traverse tree from root to tips
        # for node in self.treenode.traverse():

        #     # if children add div events
        #     if node.children:
        #         dest = min([i._schild for i in node.children])
        #         source = max([i._schild for i in node.children])
        #         time = node.height  # int(node.height)
        #         demog.add(ms.MassMigration(time, source, dest))
        #         print(dest, source, time)
        #         #^This gives me the values I need but I can't figure out 
        #         #how to get them as lists or pd.DataFrame...


    # def write_demog(self):
        # pass




#This is the .slim script structure that needs to be written by the submodule:
    #Plan is to reference pandas DF with columns: 
	#source, dest, nodeheight, gen
	# sorted by youngest --> oldest gen (or longest --> shortest nodeheight)

"""
#write beginning row:
1 { sim.addSubpop("p1", K); 
}
#for row in pandas DF, write the following:
{gen}{
    sim.addSubpopSplit("p{dest1}", K, p{source});
    sim.addSubpopSplit("p{dest2}", K, p{source});    
    p{source}.setMigrationRates(p{dest1}, 1.0);
    p{source}.setMigrationRates(p{dest2}, 1.0);    
}
{gen+1}{
    p1.setMigrationRates({dest1}, 0.0);
    p1.setMigrationRates({dest2}, 0.0);
}
#append until out of rows

#then, write last line:
self.gentime late() { 
    sim.outputFull(); 
}
"""


if __name__ == "__main__":

    import toytree
    import numpy as np

    # make a random tree with 10 tips and root height 1M
    tree = toytree.rtree.unittree(ntips=10, treeheight=1e6, seed=123)

    # set Ne values on the nodes of hte tree
    tree = tree.set_node_values(
        feature="Ne", 
        values={i: np.random.randint(10000, 100000) for i in tree.idx_dict},
    )

    # TODO: dem loads a tree and parses Ne from nodes
    dem = Demography(tree)

    # supported already
    dem = Demography(tree, Ne=10000)
    print(dem._get_demog())