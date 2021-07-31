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

#optional imports
try:
    import IPython
except ImportError:
    pass

#internal imports
from shadie.chromosome.build import ChromosomeBase


class PostSim:
    "Merges tree sequeces and recapitates"

    def __init__(
        self,
        files:list = None, #takes exactly two files (for now)
        simlength:int = None, #length in generations
        popsize:int = None, #size of each ending population
        recomb: float = None, #recomcbination rate
        mutrate:float = None, #mutations rate
        chromosome = None, #'shadie.chromosome.ChromosomeBase'
        ):
    
        """
        Reads in two SLiM .trees files, merges them, recapitates, 
        overlays neutral mutations and saves info.
        """

        self.chromosome = chromosome

        self.species = []
        ids = []
        species = []

        #read in all thre tree sequences
        for i in range(0,len(files)):
            ts = pyslim.load(files[i])
            species.append(ts)

        #merge the tree sequences
        self.merged_ts = pyslim.SlimTreeSequence(
            species[0].union(
            species[1], 
            node_mapping=[tskit.NULL for i in range(species[1].num_nodes)],
            add_populations=True,
            )
        )

        ## merge all individuals into population 1
        demographic_events = []

        for i in range(1, self.merged_ts.num_populations):
            demographic_events.append(msprime.MassMigration(
                time = simlength, source = i, destination = 0,
                proportion = 1.0))

        # OLD MSPRIME METHOD
        pop_configs = [msprime.PopulationConfiguration(initial_size=popsize)
            for _ in range(self.merged_ts.num_populations)]

        matrix = np.zeros((self.merged_ts.num_populations,self.merged_ts.num_populations))

        self.rts = self.merged_ts.recapitate(
            population_configurations=pop_configs,
            demographic_events = demographic_events,
            migration_matrix= matrix,
            recombination_rate=recomb,
            random_seed=4,
        )

        #mutate the tree sequencce
        self.mts = pyslim.SlimTreeSequence(msprime.mutate(self.rts, rate=mutrate, keep=True))

        #save pops
        alive = self.merged_ts.individuals_alive_at(0)
        num_alive = [0 for _ in range(self.merged_ts.num_populations)]
        for i in alive:
           ind = self.merged_ts.individual(i)
           num_alive[ind.population] += 1

        self.num_alive = num_alive

        edge_ids = []
        for i in range(self.mts.num_populations):
            if num_alive[i] != 0:
                edge_ids.append(i)
            else:
                pass

        self.edge_ids = edge_ids

        pop1 = []
        pop2 = []
        inds = self.mts.individuals()
        for i in range(1, inds.length):
            if inds[i].population == edge_ids[0]:
                pop1.append(inds[i].id)
            elif inds[i].population == edge_ids[1]:
                pop2.append(inds[i].id)
        self.pop1 = pop1
        self.pop2 = pop2


        #save mutation positions and which population they occurred in
        positions = []
        population = []
        for mut in self.mts.mutations():
            positions.append(int(mut.position))
            population.append(self.mts.node(mut.node).population)

        self.positions = positions
        self.population = population
        

    def simplify(
        self, 
        samplesize:int = 10, 
        random_seed = 3):
        """
        Simplify merged tree sequence. PostSim.sts and PostSim.sets
        can be accessed for tskit summary statistic functions, e.g:
        PostSim.sts.divergence(PostSim.sets)
        .diversity()
        """

        np.random.seed(random_seed)
        p1_keep_indivs = np.random.choice(self.pop1, samplesize, replace=False)
        p2_keep_indivs = np.random.choice(self.pop2, samplesize, replace=False)
        keep_indivs = np.concatenate((p2_keep_indivs, p1_keep_indivs), axis=None)
        keep_nodes = []

        for i in keep_indivs:
           keep_nodes.extend(self.mts.individual(i).nodes)
        self.sts = self.mts.simplify(keep_nodes, keep_input_roots=True)

        #save the sets
        set1 = []
        set2 = []
        for i in self.sts.individuals_alive_at(0):
            if self.sts.individual(i).population == 0:
                set1.append(i)
            else:
                set2.append(i)
        self.sets = [set1, set2]

        print(f"Before, there were {self.mts.num_samples} sample nodes "
            f"(and {self.mts.num_individuals} individuals) "
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

    
    def summary(self):
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


    def plot(self):
        """
        Plots the mutations over the chromosome as interactive altair 
        plot - but needs to be opened in vega editor (??)
        """

        if self.chromosome == None:
            logger.warning("Please init Postsim object with chromosome")
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
                    'population': self.population
                    })

                mut_positions = alt.Chart(mut_pos).mark_rule().encode(
                    color=alt.Color('population:N', 
                        scale=alt.Scale(domain=self.edge_ids, 
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
        for index, row in self.genemap.iterrows():
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



if __name__ == "__main__":
    test = PostSim(
        ["./ts0.trees", "./ts1.trees"],
        simlength = 1000,
        popsize = 1000,
        recomb = 1e-8, 
        mutrate = 1e-8)
    test.individuals()
    test.simplify()
    test.summary()
    test.plot()


