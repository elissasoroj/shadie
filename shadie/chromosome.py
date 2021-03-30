#!/usr/bin/env python

"""
Review chromosome object that will be used for a simulation. Reads in
object from subclass BuildChromosome. Otherwise, creates a single gene.

"""

#package imports
import pandas as pd
from loguru import logger
import toyplot
import altair as alt

#internal imports
from shadie.buildchromosome import Build
from shadie.mutations import MutationList
from shadie.elements import ElementList

from shadie.globals import BEN
from shadie.globals import SYN
from shadie.globals import DEL
from shadie.globals import EXON

#optional imports
try:
    import IPython
except ImportError:
    pass

#


class Chromosome:
    """
    shadie.Chromosome object stores the genome information for a SLiM
    simulation. Allows the use to inspect the genome parameters 
    before proceeding

    """

    def __init__(
        self,
        genome_size=2e3,        #will be used to calculate chromosome end (length -1)
        genome = None           #optional BuildChromosome object
        ):

        self.genome = genome
        self.gensize = genome_size
        """
        Accepts a subclass Build object to define chromosome structure.
        Otherwise, takes genome_size (length) and mutation rate and 
        creates a single gene for the simulation using default EXON 
        genomic element type and default mutation types assigned to 
        EXON (see globals.py). 

        Parameters:
        -----------

        mutation_rate(int): default = 1e-7
            Chance of mutation at each bp


        genome_size(int): default = =1e6
            Length of chromosome
        """

        if self.genome != None:
            if isinstance(self.genome, Build):
                self.mutationlist = self.genome.mutationlist
                self.elementlist = self.genome.elementlist
                self.genome = self.genome.genelements
                self.gensize = self.genome.gensize
            

        elif self.genome == None:
            g1 = [{'name': "exon", 'start': 1, 
            'finish': self.gensize - 1, 'eltype':EXON.name, 'script':EXON}]
            gene = pd.DataFrame(g1)
            self.genome = gene
            self.gensize = genome_size
            self.mutationlist = MutationList(SYN, DEL, BEN)
            self.elementlist = ElementList(EXON)
        

    def review(self, item = None):
        """
        allows user to inspect chromosome settings, including:
        mutation rates, genomic element types, and genomic elements
        also creates interactive plot of the chromosome, so the user can inspect it
        """

        if item == "mutations":
            print("Mutation Types:\n", self.mutationlist)
            #print("Mutations:\n" self.mutdict)
        elif item == "eltypes":
            print("Genomic Element Types:\n", self.elementlist)
            #print("Genomic Element Types:\n" self.eldict)
        elif item == "elements":
            print("Genomic elements:\n", self.genome)
        elif item == "chromosome":
            #Make the `rectangles` dataframe for plotting
            eltype = []
            startbase = []
            endbase = []
            y1 = []
            y2 = []
            for index, row in self.genome.iterrows():
                eltype.append(row['name'])
                startbase.append(row['start'])
                endbase.append(row['finish'])
                y1.append(0)
                y2.append(1)

            chromcoords = list(zip(eltype, startbase, endbase, y1, y2))
            rectangles = pd.DataFrame(chromcoords, columns = ['Element Type', 'x1', 'x2', 'y1', 'y2'])

            #define colors for each element
            #add will need to re-write for custom genetic elements
            color = []
            for index, row in rectangles.iterrows():
                if row["Element Type"] == "noncoding":
                    color.append("firebrick")
                elif row["Element Type"] == "exon":
                    color.append("steelblue")
                elif row["Element Type"] == "intron":
                    color.append("orange")

            #assign colors to rectangles dataframe
            rectangles.insert(5, "color", color)

            #Calculate Stats
            genecount = max(1, (rectangles["Element Type"]=='noncoding').sum()-1)
            exoncount = (rectangles["Element Type"]=='exon').sum()
            introncount = (rectangles["Element Type"]=='intron').sum()
            totexonlength = 0
            totintronlength = 0
            for index, row in rectangles.iterrows():
                if row["Element Type"]=='exon':
                    totexonlength += (1+row["x2"]-row["x1"]) 
                elif row["Element Type"]=='intron':
                    totintronlength += (1+row["x2"]-row["x1"])


            print(f"# of Genes: {genecount}\n"
                f"Average # exons per gene: {exoncount/genecount}\n"
                f"Average exon length: {totexonlength/exoncount} nt\n"
                f"Average # introns per gene: {introncount/genecount}\n"
                f"Average introns length: {totintronlength/introncount} nt\n"
                )

            print(f"Chromosome Plot:\n")
            # make the canvas and axes
            canvas = toyplot.Canvas(width=3000, height=200)
            axes = canvas.cartesian()
            axes.show = True

            #draw the rectangles
            for index, row in rectangles.iterrows():
                axes.rectangle(
                    row['x1'], row['x2'], row['y1'], row['y2'],
                    color = row['color']
                )

        elif item == "interactive":
            print("Interactive altair chromosome map:")
            #Make the `rectangles` dataframe for plotting
            eltype = []
            startbase = []
            endbase = []
            y1 = []
            y2 = []
            length = []
            for index, row in self.genome.iterrows():
                eltype.append(row['name'])
                startbase.append(row['start'])
                endbase.append(row['finish'])
                y1.append(0)
                y2.append(1)
                length.append(row['finish']-row['start'])

            chromcoords = list(zip(eltype, startbase, endbase, y1, y2, length))
            rectangles = pd.DataFrame(chromcoords, 
                columns = ['Element Type', 'x1', 'x2', 'y1', 'y2', 'length'])

            #make the altair plot
            dom = ['exon', 'intron', 'noncoding'] 
            rng = ['mediumvioletred', 'lightgoldenrodyellow', 'cornflowerblue'] 

            brush = alt.selection_interval(
                encodings=['x'], 
                mark=alt.BrushConfig(fill='red', fillOpacity = 0.700))

            ichrom = alt.Chart(rectangles).mark_rect().encode(
                x=alt.X('x1:Q', axis=alt.Axis(title='Base Pairs')),
                x2='x2:Q',
                y = alt.Y('y1:Q', axis=None),
                y2='y2:Q', 
                color=alt.Color('Element Type:N', 
                                scale=alt.Scale(domain=dom, range=rng)),
                tooltip=[
                    alt.Tooltip('Element Type', title='Element Type'),
                    alt.Tooltip('x1', title='Start'),
                    alt.Tooltip('x2', title='Stop'),
                    alt.Tooltip('length', title='Length'),
                        ]
            ).properties(
                height = 80,
                width = 850
            ).add_selection(brush)

            zoom = alt.vconcat(
                ichrom.encode(
                alt.X('x1:Q', title=None, scale=alt.Scale(domain=brush))),
                ichrom.add_selection(brush).properties(height=40),
                data=rectangles)

            print(rectangles)
            #ichrom.save('ichrom.html')
            IPython.display.display_html(zoom)
            
        else:
            logger.info("Please enter a valid item to review. Options include:\n"
                "'mutations'\n'eltypes'\n'elements'\n''chromosome'\n'interractive'")


if __name__ == "__main__":

    # generate random chromosome
    init_chromosome = Chromosome(genome_size = 2000)
    Chromosome.review(init_chromosome, item = "elements")
    Chromosome.review(init_chromosome, item = "chromosome")
    Chromosome.review(init_chromosome, item = "interactive")

    random_chromosome = Build()
    Build.random(random_chromosome)
    final_chromosome = Chromosome(genome = random_chromosome)
    final_chromosome.review("elements")
    final_chromosome.review("mutations")
