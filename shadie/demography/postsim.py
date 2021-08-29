#!/usr/bin/env python

"""
Reads in two SLiM .trees files of sister species for recapitation and 
overlaying neutral mutations.

"""

#imports
import numpy as np
import pandas as pd
import pyslim
import tskit
import msprime
import altair as alt
from loguru import logger
from typing import Union
import toytree
import toyplot
from IPython.utils import io

#optional imports
try:
    import IPython
    from IPython.display import display, SVG
except ImportError:
    pass

#internal imports
from shadie.chromosome.src.base_class import ChromosomeBase


class PostSim:
    "Merges tree sequeces and recapitates"

    def __init__(
        self,
        files:list = None, #takes exactly two files (for now)
        simlength:int = None, #length in generations
        popsize:int = None, #initial size of each population
        recomb: float = None, #recomcbination rate
        mutrate:float = None, #mutations rate
        chromosome = None, #'shadie.chromosome.ChromosomeBase'
        ):
    
        """
        Reads in two SLiM .trees files, merges them, recapitates, 
        overlays neutral mutations and saves info.
        """

        self.chromosome = chromosome
        self.simlength = simlength
        self.mutrate = mutrate/2
        self.recomb = recomb
        self.popsize = popsize
        self.pops = None

        self.species = []
        ids = []
        species = []

        #read in all thre tree sequences
        for i in range(0,len(files)):
            ts = pyslim.load(files[i])
            species.append(ts)

        #merge the p0 and p1 populations in both edges
        tablelist = []
        mod_tslist = []
        for ts in species:
            tables = ts.tables
            tables.nodes.population = np.zeros(tables.nodes.num_rows, dtype=np.int32)
            modts = pyslim.load_tables(tables)
            mod_tslist.append(modts)

        #remove extra population
        onepop_tslist = []
        for ts in mod_tslist:
            onepop = ts.simplify(keep_input_roots=True,
                keep_unary_in_individuals = True)
            onepop_tslist.append(onepop)
        
        self.onepop_tslist = onepop_tslist

    def fullworkflow(self):

        if len(self.onepop_tslist) > 1:
            self.merge(self.onepop_tslist[0], self.onepop_tslist[1])
            self.recapitate(self.merged_ts)
            self.mutate(self.rts)
        else:
            self.recapitate(self.onepop_tslist[0])
            self.mutate(self.rts)

    
    def merge(self, ts1, ts2):    #merge the tree sequences
        self.merged_ts = pyslim.SlimTreeSequence(
            ts1.union(
            ts2, 
            node_mapping=[tskit.NULL for i in range(ts2.num_nodes)],
            add_populations=True,
            )
        )

        #save pops
        alive = self.merged_ts.individuals_alive_at(0)
        num_alive = [0 for _ in range(self.merged_ts.num_populations)]
        for i in alive:
           ind = self.merged_ts.individual(i)
           num_alive[ind.population] += 1

        self.num_alive = num_alive

        edge_ids = []
        for i in range(self.merged_ts.num_populations):
            if num_alive[i] != 0:
                edge_ids.append(i)
            else:
                pass

        self.edge_ids = edge_ids

        pop1 = []
        pop2 = []
        inds = self.merged_ts.individuals()
        for i in range(1, inds.length):
            if inds[i].population == self.edge_ids[0]:
                pop1.append(inds[i].id)
            elif inds[i].population == self.edge_ids[1]:
                pop2.append(inds[i].id)
        
        self.pops = [pop1, pop2]

        #save starting pops
        # startalive = self.merged_ts.individuals_alive_at(self.simlength-1)
        # num_startalive = [0 for _ in range(self.merged_ts.num_populations)]
        # for i in alive:
        #    ind = self.merged_ts.individual(i)
        #    num_alive[ind.population] += 1

        # self.num_startalive = num_startalive


    def recapitate(self, ts):    ## merge all individuals into population 1
        """
        recapitate the treesequence. This function drops any empty 
        nodes due to haploid genomes 
        """
        
        #drop empty nodes
        nodeslist = []
        for i in ts.edges():
            nodeslist.append(ts.edge(i.id).parent)
            nodeslist.append(ts.edge(i.id).child)

        keep_nodes = []
        [keep_nodes.append(x) for x in nodeslist if x not in keep_nodes]
        len(keep_nodes)

        ts = ts.simplify(samples = keep_nodes, keep_input_roots=True)

        demographic_events = []
        if ts.num_populations > 1:
            for i in range(1, ts.num_populations):
                demographic_events.append(msprime.MassMigration(
                    time = self.simlength+10, source = i, destination = 0,
                    proportion = 1.0))

        # OLD MSPRIME METHOD
        # pop_configs = []
        # for pop in ts.populations():
        #     pop_configs.append(msprime.PopulationConfiguration(
        #         initial_size = 1+self.num_startalive[pop.id]))
        
        pop_configs = [msprime.PopulationConfiguration(
            initial_size=self.popsize)
            for _ in range(ts.num_populations)]

        matrix = np.zeros((ts.num_populations, ts.num_populations))

        self.rts = ts.recapitate(
            population_configurations=pop_configs,
            demographic_events = demographic_events,
            migration_matrix= matrix,
            recombination_rate=self.recomb,
            random_seed=4,
        )

    def mutate(self, ts):    
        "mutate the tree sequence"
        
        self.mts = pyslim.SlimTreeSequence(msprime.mutate(
            ts, rate=self.mutrate, keep=True))

        #save mutation positions and which population they occurred in
        positions = []
        popids = []
        allpositions = []
        for mut in self.mts.mutations():
            allpositions.append(int(mut.site))
            popids.append(self.mts.node(mut.node).population)
            if mut.derived_state != '1':
                positions.append(int(mut.site))

        self.positions = positions
        self.allpositions = allpositions
        self.popids = popids

    def simplify(
        self, 
        ts=None,
        samplesize:int = 10, 
        random_seed:Union[int, None]=None):
        """
        Simplify merged tree sequence. Note that this simplified tree
        sequence cannot be recapitated.
        PostSim.sts and PostSim.sets can be accessed for tskit summary 
        statistic functions, e.g:
        PostSim.sts.divergence(PostSim.sets)
        .diversity()
        """
        keep_indivs = []
        np.random.seed(random_seed)
       
        if self.pops is not None:
            ts = self.mts
            for pop in self.pops:
                inds = np.random.choice(pop, samplesize, replace=False)
                for i in inds:
                    keep_indivs.append(i)

        else:
            inds = []
            for i in ts.individuals():
                inds.append(i.id)
            keep_indivs = np.random.choice(inds, samplesize, replace=False)

        keep_nodes = []
        for i in keep_indivs:
           keep_nodes.extend(ts.individual(i).nodes)
        
        self.sts = ts.simplify(keep_nodes)

        #save the sets
        set1 = []
        set2 = []
        for i in self.sts.individuals_alive_at(0):
            if self.sts.individual(i).population == 0:
                set1.append(i)
            else:
                set2.append(i)
        self.sets = [set1, set2]

        print(f"Before, there were {ts.num_samples} sample nodes "
            f"(and {ts.num_individuals} individuals) "
            f"in the tree sequence, and now there are {self.sts.num_samples} "
            f"(sample nodes and {self.sts.num_individuals} individuals).\n\n"
            "Use `.sts` to access simplified tree sequnce and  `.sets`"
            "for summary statistics")
    
    
    def individuals(self):
        """
        Returns number of individuals in each population
        """

        for pop, num in enumerate(self.num_alive):
           print(f"Number of individuals in population {pop}: {num}")

    
    def quickplot(self, ts, node, size:tuple = (800,400)):
        display(SVG(ts.at_index(node).draw_svg(size = size, y_axis=True,)))

    def treeplot(self, ts, node, width:int = 600, height:int = 700):
        ts_newick = ts.at_index(node).newick()

        modtree = toytree.tree(ts_newick)
        canvas, axes, mark = modtree.draw(width=width, height=height)

        # show the axes coordinates
        axes.show = True
        axes.x.ticks.show = True
        axes.y.ticks.show = True

        axes.vlines(1-self.simlength, style={"stroke": "blue"});


    def stats(
        self, 
        ts=None, 
        samplesets:list = None, 
        samplesize:int = 20,
        sampletimes:int = 50):
        """
        View a summary of tskit statistics on the treesequence. 
        If a tree sequence and samples sets are provided, stats() will
        calculate summary statistics using sample set.
        If only a tree sequence is provided, the tree sequence will be
        sampled sampletimes times (Default = 50) at samplesize (Default = 20).
        If no tree sequence and sample sets are provided, stats() will 
        attempt to use simplified PostSim.mts sequence. 
        """

        if samplesets is not None:
           
            print(
            f"Divergence: {ts.divergence(samplesets)}\n"
            f"Diversity: {ts.diversity(samplesets)}\n"
            f"Fst: {ts.Fst(samplesets)}\n"
            f"Tajima's D: {ts.Tajimas_D(samplesets)}\n"
            f"Root age: {ts.max_root_time}"
            )

        else:
            if ts is not None:
                pass
            elif self.mts is not None:
                ts = self.mts 
            else:
                logger.warning("Please input a tree sequence")

            divg = []
            divs1 = []
            divs2 = []
            fst = []
            tajd1 = []
            tajd2 = []

            with io.capture_output() as captured:
                for i in range (0, sampletimes):
                    self.simplify(ts, samplesize)
                    divg.append(self.sts.divergence(self.sets))
                    divs = self.sts.diversity(self.sets)
                    if isinstance(divs[0], float):
                        divs1.append(divs[0]) 
                    if isinstance(divs[1], float):
                        divs2.append(divs[1]) 
                    fst.append(self.sts.Fst(self.sets))
                    tajds = self.sts.Tajimas_D(self.sets)
                    if isinstance(tajds[0], float):
                        tajd1.append(tajds[0]) 
                    if isinstance(tajds[1], float):
                        tajd2.append(tajds[1]) 

            print(
                "Sampled 20 individuals from each population "
                f"{sampletimes} times\n"
                f"Divergence: {np.mean(divg)}\n"
                f"Diversity: pop1 = {np.mean(divs1)}, pop2 = {np.mean(divs2)}\n"
                f"Fst: {np.mean(fst)}\n"
                f"Tajima's D: pop1 = {np.mean(tajd1)}, pop2 = {np.mean(tajd2)}\n"
                f"Root age: {ts.max_root_time}"
                )


    def mutsummary(self):
        "View selection coefficient distributions and mutation type counts"

        #calculate number of neutral mutations
        neutral = 0
        for mut in self.mts.mutations():
            if mut.derived_state == '1':
                neutral += 1

        #record fitness coefficients
        coeffs = []
        for mut in self.mts.mutations():
            mut_list = mut.metadata["mutation_list"]
            for i in mut_list:
                coeffs.append(i["selection_coeff"])
        
        self.coeffs = coeffs

        #altair plot of mutation distrbutions
        mean = sum(coeffs)/len(coeffs)
        dist = np.histogram(coeffs, 50)
        source = pd.DataFrame({'counts': dist[0], 'values': dist[1][:-1], 
            'mean':mean}, columns=['values', 'counts', 'mean']) 
        base = alt.Chart(source)

        histo = base.mark_bar(
            opacity=0.5,
            color = alt.ColorName("cornflowerblue")
        ).encode(
            alt.X('values:Q', axis=alt.Axis(title='Fitness Effect')),
            alt.Y('counts:Q'),
            color=alt.condition(
                alt.datum.values > 0,
                alt.value("lightseagreen"),  # The positive color
                alt.value("indianred")  # The negative color
            ),
            tooltip = [
                alt.Tooltip('values', title='Fitness Effect'),
            ]
        )
        mean = base.mark_rule(color = alt.ColorName("goldenrod"),
                              opacity=0.4).encode(
            x='mean:Q',
            size=alt.value(3),
            tooltip=[
                alt.Tooltip('mean(values)', title='Mean Fitness Effect'),
            ])

        print(
            "Note: this plot does not contain neutral mutations overlaid "
            "with `msprime`."
            )
        IPython.display.display_html(histo + mean)
        print(
            f"Total mutations: {self.mts.num_mutations}\n"
            f"Neutral mutations: {neutral}\n"
            f"Non-neutral mutations: {self.mts.num_mutations - neutral}\n"
            )


    def plotmuts(self, chromosome=None):
        """
        Plots the mutations over the chromosome as interactive altair 
        plot - but needs to be opened in vega editor (??)
        """

        if self.chromosome == None:
            if chromosome is not None:
                self.chromosome = chromosome
            else:
                logger.warning("PostSim object was not initialised with"
                    "a chromosome object - please provide one for plot"
                    "function.")
        else:
            self.chromosome.inspect()
            ichrom = self.chromosome.ichrom

            alt.data_transformers.disable_max_rows()
            brush = alt.selection_interval(
                    encodings=['x'], 
                    mark=alt.BrushConfig(fill='red', fillOpacity = 0.700))

            fadedchrom = ichrom.mark_rect(opacity = 0.4)


            if len(self.edge_ids) > 1:
                mut_pos = pd.DataFrame({
                    'x': self.positions,
                    'population': str("pop" + self.popids)
                    })

                mut_positions = alt.Chart(mut_pos).mark_rule().encode(
                    color=alt.Color('population', 
                        scale=alt.Scale(domain= self.edge_ids, 
                            range=["mediumblue", "crimson"])),
                    x='x:Q',
                    size=alt.value(1),
                    tooltip=[alt.Tooltip('x', title='Position'),]
                    )
            else:
                mut_pos = pd.DataFrame({'x': self.positions})
                self.mut_pos = mut_pos
                mut_positions = alt.Chart(mut_pos).mark_rule(color = alt.ColorName("mediumblue")).encode(
                        x='x:Q',
                        size=alt.value(1),
                        tooltip=[
                            alt.Tooltip('x', title='Position'),
                        ])

            layered = alt.layer(fadedchrom, mut_positions)

            zoomtest = alt.vconcat(
                layered.encode(
                alt.X('x1:Q', title=None, scale=alt.Scale(domain=brush))).properties(height = 80),
                layered.add_selection(brush).properties(height=30))

            print("Note: this plot must be opened in the Vega editor for "
                "interactive features to work")
            IPython.display.display_html(zoomtest)

            self.plot = zoomtest



    def dNdS(self):
        "calculate dN/dS"
        ranges = []
        for index, row in self.chromosome.data.iterrows():
            ranges.append(range(row['exonstart'], row['exonstop']))
            
        c_syn = []
        c_non = []
        nc_neut = []
        for mut in self.tscoal.mutations():
            if mut.derived_state != '1':
                c_non.append(mut.id)      
            elif mut.derived_state == '1':
                position = int(mut.position)
                for r in ranges:
                    if position in r:
                        c_syn.append(mut.id)
                    elif position not in r:
                        nc_neut.append(mut.id)
        
        dnds = len(c_non)/len(c_syn)
        print(f"N:S for the coding regions in whole genome = {dnds}")


    def mktest(self):
        "Perform mk-test"
        pass

class NeutralSim:
    "tools for a SLiM simulation in which neutral mutations are included."

    def __init__(
        self,
        files:list = None, #only two files for now
        simlength:int = None, #length in generations
        #popsize:int = None, #size of each ending population
        recomb: float = None, #recomcbination rate
        mutrate:float = None, #mutations rate
        chromosome = None, #'shadie.chromosome.ChromosomeBase'
        ):

        self.chromosome = chromosome
        self.files = files
        self.recomb = recomb
        self.simlength = simlength
        self.mutrate = mutrate/2

    def read(self):
        self.ts = pyslim.load(self.files[i])

    def merge(self):
        "Merges tree sequences"

        self.species = []
        ids = []
        species = []

        #read in all thre tree sequences
        for i in range(0,len(self.files)):
            ts = pyslim.load(self.files[i])
            species.append(ts)

        self.merged_ts = pyslim.SlimTreeSequence(
            species[0].union(
            species[1], 
            node_mapping=[tskit.NULL for i in range(species[1].num_nodes)],
            add_populations=True,
            )
        )


    def mutate(self):
        "mutates the treesequence"


    def recapitate(self, treeseq):
        self.treeseq = treeseq
        #save pops

        alive = self.treeseq.individuals_alive_at(0)
        num_alive = [0 for _ in range(self.treeseq.num_populations)]
        for i in alive:
           ind = self.treeseq.individual(i)
           num_alive[ind.population] += 1

        self.num_alive = num_alive

        edge_ids = []
        for i in range(self.treeseq.num_populations):
            if num_alive[i] != 0:
                edge_ids.append(i)
            else:
                pass

        # OLD MSPRIME METHOD
        pop_configs = []
        for pop in self.treeseq.populations():
            pop_configs.append(msprime.PopulationConfiguration(
                initial_size = 1+self.num_alive[pop.id]))
        
        # pop_configs = [msprime.PopulationConfiguration(initial_size=popsize)
        #     for _ in range(self.merged_ts.num_populations)]

        matrix = np.zeros((self.treeseq.num_populations,self.treeseq.num_populations))

        rts = self.treeseq.recapitate(
            population_configurations=pop_configs,
            migration_matrix= matrix,
            recombination_rate=self.recomb,
            random_seed=4,
        )

        self.rts = rts

    def mutcount(self):
        "Returns counts of neutral and non-neutral mutations"

        nonneutral = []
        neutral = []
        mutationcount = 0
        neutralcount = 0

        for mut in self.mts.mutations():
            oldest = len(mut.metadata["mutation_list"]) - 1
            if mut.metadata["mutation_list"][oldest]["selection_coeff"] != 0.0:
                nonneutral.append(mut)
                mutationcount += 1
            else:
                neutral.append(mut)
                neutralcount += 1

        self.neutral = neutral
        self.nonneutral =  nonneutral
        self.mutationcount = mutationcount
        self.neutralcount = neutralcount

        print(f"The tree sequence has {self.neutralcount} neutral mutations and "
            f"{self.mutationcount} non-neutral mutations")



if __name__ == "__main__":
    test = PostSim(
        ["./ts0.trees", "./ts1.trees"],
        simlength = 1000,
        popsize = 1000,
        recomb = 1e-8, 
        mutrate = 1e-8)
    test.fullworkflow()
    test.stats()

