#!/usr/bin/env python

"""
Performs post-simulation calculations
"""

#imports
import numpy as np
import pandas as pd
import pyslim
import altair as alt
import toyplot

#optional imports
try:
    import IPython
except ImportError:
    pass

#internal imports
from shadie.chromosome import Chromosome
from shadie.main import Shadie
from shadie.slimcoal import Coal


class PostSim:
    """
    Post-simulation summary statistics
    """
    def __init__(
        self,
        shadie,
        ):
    
        """
        Builds script to run SLiM3 simulation

        Parameters:
        -----------
        shadie: (Shadie class object)
            Reads in Shadie class object
        """
        if isinstance(shadie, Shadie):
            self.shadie = shadie

        self.tsraw = pyslim.load(self.shadie.outname)

        tsout = Coal(self.tsraw)
        tsout.slimcoal()

        self.tscoal = tsout.tscoal
        self.genemap = shadie.genemap
        self.genome = shadie.genome

        #positions
        positions = []
        for mut in self.tscoal.mutations():
            positions.append(int(mut.position))
        self.positions = positions

    
    def summary(self):
        "view selection coefficient distributions and mutation type counts"
        #calculate number of neutral mutations
        neutral = 0
        for mut in self.tscoal.mutations():
            if mut.derived_state == '1':
                neutral += 1

        #record fitness coefficients
        coeffs = []
        for mut in self.tscoal.mutations():
            mut_list = mut.metadata["mutation_list"]
            for i in mut_list:
                coeffs.append(i["selection_coeff"])
        
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
            f"Total mutations: {self.tscoal.num_mutations}\n"
            f"Neutral mutations: {neutral}\n"
            f"Non-neutral mutations: {self.tscoal.num_mutations - neutral}\n"
            )

        #static toyplot
        print("Mutation positions along chromosome:")
        chrom = Chromosome(self.genome)
        chrom.toyplot()
        self.rectangles = chrom.rectangles

        canvas = toyplot.Canvas(width=2500, height=200)
        axes = canvas.cartesian()
        axes.show = False

        #draw the rectangles
        for index, row in self.rectangles.iterrows():
            axes.rectangle(
                row['x1'], row['x2'], row['y1'], row['y2'],
                color = row['color'],
                style={"opacity":0.6},
            )
        #draw the positions
        lines = axes.vlines(self.positions, style={"stroke":"blue", "stroke-width":2})


    def zoomplot(self):
        "interactive altair plot - but needs to be opened in vega editor"
        chrom = Chromosome(self.genome)
        chrom.altair()
        ichrom = chrom.ichrom

        brush = alt.selection_interval(
                encodings=['x'], 
                mark=alt.BrushConfig(fill='red', fillOpacity = 0.700))

        fadedchrom = ichrom.mark_rect(opacity = 0.4)

        mut_pos = pd.DataFrame({'x': self.positions})
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
        print(f"dN/dS for the whole genome = {dnds}")


    def mktest(self):
        "Perform mk-test"
        pass

    def pointmutations(self):
        "lists point mutations (for nucleotide simulation only)"
        self.ts = self.tscoal
        neutral = 0
        M = [[0 for _ in pyslim.NUCLEOTIDES] for _ in pyslim.NUCLEOTIDES]
        for mut in self.ts.mutations():
            mut_list = mut.metadata["mutation_list"]
            if mut_list == []:
                neutral += 1
            else:
                k = np.argmax([u["slim_time"] for u in mut_list])
                derived_nuc = mut_list[k]["nucleotide"]
                if mut.parent == -1:
                    acgt = self.tsraw.reference_sequence[int(mut.position)]
                    parent_nuc = pyslim.NUCLEOTIDES.index(acgt)
                else:
                    parent_mut = self.ts.mutation(mut.parent)
                    assert(parent_mut.site == mut.site)
                    parent_nuc = parent_mut.metadata["mutation_list"][0]["nucleotide"]
                M[parent_nuc][derived_nuc] += 1
                
        print("{}\t{}\t{}".format('ancestr', 'derived', 'count'))
        for j, a in enumerate(pyslim.NUCLEOTIDES):
            for k, b in enumerate(pyslim.NUCLEOTIDES):
                print("{}\t{}\t{}".format(a, b, M[j][k]))
        print(f"\nneutral mutations overlaid with coalescent: {neutral}")

    def aminoacids (self):
        "lists amino acid mutations (for nucleotide simulation only)"
        self.ts = self.tscoal

        slim_gen = self.ts.metadata["SLiM"]["generation"]
        
        M = np.zeros((4,4,4,4), dtype='int')
        for mut in self.ts.mutations():
            pos = self.ts.site(mut.site).position
            # skip mutations at the end of the sequence
            if pos > 0 and pos < self.ts.sequence_length - 1:
                mut_list = mut.metadata["mutation_list"]   
                k = np.argmax([u["slim_time"] for u in mut_list])
                derived_nuc = mut_list[k]["nucleotide"]
                left_nuc = self.ts.nucleotide_at(mut.node, pos - 1,
                time = slim_gen - mut_list[k]["slim_time"] - 1.0)
                right_nuc = self.ts.nucleotide_at(mut.node, pos + 1,
                time = slim_gen - mut_list[k]["slim_time"] - 1.0)
                if mut.parent == -1:
                    acgt = self.ts.reference_sequence[int(mut.position)]
                parent_nuc = pyslim.NUCLEOTIDES.index(acgt)
            else:
                parent_mut = self.ts.mutation(mut.parent)
                assert(parent_mut.site == mut.site)
                parent_nuc = parent_mut.metadata["mutation_list"][0]["nucleotide"]
        
            M[left_nuc, parent_nuc, right_nuc, derived_nuc] += 1
        print("{}\t{}\t{}".format('ancestral', 'derived', 'count'))
        for j0, a0 in enumerate(pyslim.NUCLEOTIDES):
            for j1, a1 in enumerate(pyslim.NUCLEOTIDES):
                for j2, a2 in enumerate(pyslim.NUCLEOTIDES):
                    for k, b in enumerate(pyslim.NUCLEOTIDES):
                        print("{}{}{}\t{}{}{}\t{}".format(a0, a1, a2, a0, b, a2,
                        M[j0, j1, j2, k]))



if __name__ == "__main__":
    pass