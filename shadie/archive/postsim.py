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
    "Post-simulation analysis and summary statistics"

    def __init__(
        self,
        shadie = None,
        file = None,
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
            self.genemap = shadie.genemap
            self.genome = shadie.genome

        elif shadie is None:
            self.tsraw = pyslim.load(file)

        tsout = Coal(self.tsraw)
        tsout.slimcoal()

        self.tscoal = tsout.tscoal

        #positions
        positions = []
        for mut in self.tscoal.mutations():
            positions.append(int(mut.position))
        self.positions = positions

    
    def summary(self):
        "View selection coefficient distributions and mutation type counts"
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
        chrom = Chromosome(genome = self.genome)
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
        print(f"N:S for the coding regions in whole genome = {dnds}")


    def mktest(self):
        "Perform mk-test"
        pass



if __name__ == "__main__":
    pass
