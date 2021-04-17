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

#optional imports
try:
    import IPython
    from IPython.display import display
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
        genome_size=2e5,        #will be used to calculate chromosome end (length -1)
        genome = None           #optional BuildChromosome object
        ):
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
        self.genome = genome

        if self.genome is not None:
            if isinstance(self.genome, Build):
                self.mutationlist = genome.mutationlist
                self.elementlist = genome.elementlist
                self.genome = genome.genelements
                self.gensize = genome.gensize
            

        elif self.genome is None:
            onegene = Build(genecount = 1, genome_size = genome_size)
            onegene.genes()
            self.mutationlist = onegene.mutationlist
            self.elementlist = onegene.elementlist
            self.genome = onegene.genelements
            self.gensize = onegene.gensize
        

    def altair(self):
        "Makes interactive altair plot"
        eltype = []
        altname = []
        startbase = []
        endbase = []
        y1 = []
        y2 = []
        length = []
        for index, row in self.genome.iterrows():
            eltype.append(row['type'])
            if pd.isna(row['name']):
                altname.append(row['type'])
            else: 
                altname.append(row['name'])
            startbase.append(row['start'])
            endbase.append(row['finish'])
            y1.append(0)
            y2.append(1)
            length.append(row['finish']-row['start'])

        chromcoords = list(zip(eltype, altname, startbase, endbase, y1, y2, length))
        rectangles = pd.DataFrame(chromcoords, 
            columns = ['type', 'altname', 'x1', 'x2', 'y1', 'y2', 'length'])

        #set the colors
        extypecount = rectangles.loc[rectangles['type']=='exon'].altname.nunique()
        extypes = list(rectangles.loc[rectangles['type']=='exon'].altname.unique())
        intypecount = rectangles.loc[rectangles['type']=='intron'].altname.nunique()
        intypes = list(rectangles.loc[rectangles['type']=='intron'].altname.unique())
        nctypecount = rectangles.loc[rectangles['type']=='noncoding'].altname.nunique()
        nctypes = list(rectangles.loc[rectangles['type']=='noncoding'].altname.unique())

        ncodcolors = ['mediumvioletred', 'lightcoral','firebrick', 'crimson', 'lightpink']
        incolors = ['lemonchiffon', 'gold', 'orange',  'yellow', 'khaki']
        excolors = ['cornflowerblue', 'mediumblue','dodgerblue', 'darkslateblue', 'skyblue']


        dom = extypes + intypes + nctypes
        rng = excolors[0:extypecount] + incolors[0:intypecount] + ncodcolors[0:nctypecount]

        brush = alt.selection_interval(
            encodings=['x'], 
            mark=alt.BrushConfig(fill='red', fillOpacity = 0.700))

        #make the altair plot
        ichrom = alt.Chart(rectangles).mark_rect().encode(
            x=alt.X('x1:Q', axis=alt.Axis(title='Base Pairs')),
            x2='x2:Q',
            y = alt.Y('y1:Q', axis=None),
            y2='y2:Q', 
            color=alt.Color('altname:N', 
                            scale=alt.Scale(domain=dom, range=rng)),
            tooltip=[
                alt.Tooltip('type', title='Element Type'),
                alt.Tooltip('altname', title='Name'),
                #alt.Tooltip('mutations', title='Mutations'),
                alt.Tooltip('length', title='Length'),
                alt.Tooltip('x1', title='Start'),
                alt.Tooltip('x2', title='Stop'),
                    ]
        ).properties(
            width = 850
        ).add_selection(brush)

        zoom = alt.vconcat(
            ichrom.encode(
            alt.X('x1:Q', title=None, scale=alt.Scale(domain=brush))
            ).properties(height = 80),
            ichrom.add_selection(brush).properties(height=40),
            data=rectangles)

        #zoom.save('zoom.html')
        self.ichrom = ichrom    #for postsim
        self.zoom = zoom        #for interactive plot

    def toyplot(self):
        "Makes static toyplot"
        eltype = []
        startbase = []
        endbase = []
        y1 = []
        y2 = []
        for index, row in self.genome.iterrows():
            eltype.append(row['type'])
            startbase.append(row['start'])
            endbase.append(row['finish'])
            y1.append(0)
            y2.append(1)

        chromcoords = list(zip(eltype, startbase, endbase, y1, y2))
        rectangles = pd.DataFrame(chromcoords, columns = [
            'Element Type', 'x1', 'x2', 'y1', 'y2'])

        #define colors for each element
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
        self.genecount = max(1, (rectangles["Element Type"]=='noncoding').sum()-1)
        self.exoncount = (rectangles["Element Type"]=='exon').sum()
        self.introncount = (rectangles["Element Type"]=='intron').sum()
        self.totexonlength = 0
        self.totintronlength = 0
        for index, row in rectangles.iterrows():
            if row["Element Type"]=='exon':
                self.totexonlength += (1+row["x2"]-row["x1"]) 
            elif row["Element Type"]=='intron':
                self.totintronlength += (1+row["x2"]-row["x1"])

        if self.introncount == 0:
            self.avintron = 0
        else:
            self.avintron = self.totintronlength/self.introncount

        self.rectangles = rectangles

    def review(self, item = None):
        """
        allows user to inspect chromosome settings, including:
        mutation rates, genomic element types, and genomic elements
        also creates interactive plot of the chromosome, so the user can inspect it
        """

        if item == "mutations":
            print('\033[1m' + "Mutation Types:\n" + '\033[0m', self.mutationlist, "\n")
            #print("Mutations:\n" self.mutdict)
        elif item == "eltypes":
            print('\033[1m' + "Genomic Element Types:\n" + '\033[0m', self.elementlist, "\n")
            #print("Genomic Element Types:\n" self.eldict)
        elif item == "elements":
            df = pd.DataFrame(self.genome)
            print('\033[1m' + "Genomic Elements:\n" + '\033[0m')
            display(df)

        elif item == "chromosome":
            self.toyplot()
            
            print(
                '\033[1m' + "Chromosome Summary\n" + '\033[0m'
                f"# of Genes: {self.genecount}\n"
                f"Average # exons per gene: {self.exoncount/self.genecount}\n"
                f"Average exon length: {self.totexonlength/self.exoncount} nt\n"
                f"Average # introns per gene: {self.introncount/self.genecount}\n"
                f"Average introns length: {self.avintron} nt\n"
                )

            print(f"Static Chromosome Plot:\n")
            """
            this is a simplified plot showing exons, coding regions, 
            and non-codding regions collapsed, even if they consist of 
            more that one genomic element type
            """
            # make the canvas and axes with toyplot
            canvas = toyplot.Canvas(width=2400, height=200)
            axes = canvas.cartesian()
            axes.show = True

            #draw the rectangles
            for index, row in self.rectangles.iterrows():
                axes.rectangle(
                    row['x1'], row['x2'], row['y1'], row['y2'],
                    color = row['color']
                )

        elif item == "interactive":
            self.altair()
            print("Interactive altair chromosome map:")
            IPython.display.display_html(self.zoom)
            
        else:
            logger.info("Please enter a valid item to review. Options include:\n"
                "'mutations'\n'eltypes'\n'elements'\n''chromosome'\n'interractive'")


if __name__ == "__main__":

    # generate random chromosome
    init_chromosome = Chromosome(genome_size = 2000)
    Chromosome.review(init_chromosome, item = "elements")
    init_chromosome.review("interactive")
    init_chromosome.review("chromosome")
