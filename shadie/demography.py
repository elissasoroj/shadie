#!/usr/bin/env python

"""
Traverses phylogeny from tips to root and builds corresponding population demography for SLiM

"""

#imports
import toytree
import pandas as pd

class Demography:
    def __init__(self, tree, gens=2000):
        
        # store input params
        self.tree = tree #.tree object from toytree
        self.gentime = gens #number of generations simulation will run
        self.treeheight = #read treeheight from tree object
        
        # will be used to store output results
        self.demog = None #demography in proper format

    def _get_demog(self):
        """
        Returns demography scenario based on an input tree and admixture
        edge list with events in the format (source, dest, start, end, rate).
        Time on the tree is defined in units of generations.
        """
# Define demographic events for msprime
        demogdf = pd.DataFrame()
    

        # traverse tree from root to tips
        for node in self.treenode.traverse():
            if node.children:
                node._schild = min([i.idx for i in node.get_descendants()])
            else:
                node._schild = node.idx

        # traverse tree from root to tips
        for node in self.treenode.traverse():

            # if children add div events
            if node.children:
                dest = min([i._schild for i in node.children])
                source = max([i._schild for i in node.children])
                time = node.height  # int(node.height)
                demog.add(ms.MassMigration(time, source, dest))
                print(dest, source, time)
                #^This gives me the values I need but I can't figure out 
                #how to get them as lists or pd.DataFrame...


#This is the slim. script structure that needs to be written by the class:
#Plan is to reference pandas DF with columns: 
	#source, dest, nodeheight, gen
	# sorted by youngest --> oldest gen (or longest --> shortest nodeheight)

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

