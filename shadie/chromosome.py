#!/usr/bin/env python

"""
Generates script for SLiM simulation

"""

#imports
import math
import io
import pandas as pd
import numpy as np
from loguru import logger
import toyplot
import altair as alt

#optional imports
try:
    import IPython 
except ImportError:
    pass

#standard mutation types, with structure:
#(name, dominance, distribution, {following depends on distribution})


NEUT = '("m100", 0.5, "f", 0.0)'      #neutral mutation
SYN = '("m101", 0.5, "f", 0.0)'           #synonymous
DEL = '("m102", 0.1, "g", -0.03, 0.2)'    #deleterious
BEN = '("m103", 0.8, "e", 0.1)'           #beneficial

#standard genetic elements have the structure:
#(name, mutations, frequency of each mutation)
EXON = '("g100", c(m100,m101,m102), c(2,8,0.1))'  #exon
INTRON = '("g101", c(m100,m102), c(9,1))'         #intron
NONCOD = '("g102", c(m100), 1)'               #non-coding


class Chromosome:
    """
    Chromosome object created based on user parameters
    """

    def __init__(
        self,
        chromtype = "random",       #"random" or "dict" 
        genes=None,             #number of genes on the chromosome (if None, random)
        introns = None,         #number of introns per gene (if None, random)
        exons = None,           #number of exons per gene (if None, random)
        mutation_rate = 1e-7,   #mutation rate will be used to calculate mutation matrix
        genome_size=1e6,        #will be used to calculate chromosome end (length -1)
        ):
    
        self.type = chromtype
        self.genes = genes
        self.introns = introns
        self.exons = exons
        self.mutrate = mutation_rate
        self.gensize = genome_size
        """
        Builds the chromosome for SLiM3 simulation

        Parameters:
        -----------
        chromtype(str): default = "random"
            "random" will generate a random chromosome
            "dict" will accept a dictionary object

        genes(int): default = None
            Number of genes on chromosome
            Can be supplied as `int` for random chromosome generation
            or as dict object for explicit chromosome generation, e.g:

            dict(name = "name", mutations = (fitness effect,  dominance), start = bp, end = bp)
            dict(name = "a", mutations = (-.01, 0.5), start = 2000, end = 3000)

            # call simulate with details on genome structure (and which life stage selection occurs?)
            mod.simulate(
            dict(name='a', selection=-0.01, start=2000, end=3000, haploid=True),
            dict(name='b', selection=-0.02, start=5000, end=6000, haploid=True),
            dict(name='c', selection=-0.03, start=9000, end=10000, haploid=True),
            dict(name='A', selection=0.01, start=2000, end=3000, haploid=False),
            dict(name='B', selection=0.02, start=5000, end=6000, haploid=False),
            dict(name='C', selection=0.03, start=9000, end=10000, haploid=False),
            )

        introns (int): default = None
            For random generator: number of introns per gene; default = random
            For dict generator: number of introns per gene; default = 0

        exons (int): default = None
           For random generator: number of exons per gene; default = random
           For dict generator: number of exons per gene; default = 1

        mutation_rate(int): default = 1e-7
            Chance of mutation at each bp


        genome_size(int): default = =1e6
            Length of chromosome
        """


    def make(self, *args):
        """
        argument must be in form of dictionary objects in the format:
        a = dict(name = "a", mutations = (-.01, 0.5), start = 2000, end = 3000)
        OR
        .csv file
        OR 
        pandas DataFrame
        """
        if self.type == "custom":
            if isinstance(arg, io.TextIOBase):
                #save .csv as dataframe
                pass

            elif isinstance(args, pd.DataFrame):
                self.genedf = arg

            else:
                genelist = []
                for i in args:
                    if isinstance(i, dict):
                        genelist.append(i)
                genedf = pd.DataFrame(genelist)

        elif self.type == "random":
                    pass

    def generate(self):
        "generates chromosome based on explicit user structure"
        for gene in self.genedf:
            #make the mutation types
            mutdict = {}
            #...
            self.mutdict = mutdict

            #make the genomic element types
            elemdict = {}
            #...
            self.selemdict = elemdict

            #write the initializeGenomicElement lines to a dict
            initdict = {}
            #...
            self.initdict = initdict


    def make_rand(self):
        "generates a random chromosome"
        if self.type == "random":
            if self.exons == None and self.introns == None:
                genelements = pd.DataFrame(
                columns=['name', 'start', 'finish', 'eltype'],
                data=None,
                )
                base = int(0)
                finalnc_length = np.random.randint(100, 5000)
                end = self.gensize - finalnc_length
                logger.debug("Made objects: base = {base}, genelements = {}, end = {end}")

                while (base < end):
                    #make initial non-coding region
                    nc_length = np.random.randint(100, 5000)
                    genelements.loc[base, 'name'] = "noncoding"
                    genelements.loc[base, 'eltype'] = NONCOD
                    genelements.loc[base, 'start'] = base
                    genelements.loc[base, 'finish'] = base + nc_length - 1
                    base = base + nc_length
                    
                    #make first exon
                    ex_length = round(np.random.lognormal(math.log(250), math.log(1.3))) + 1
                    genelements.loc[base, 'name'] = "exon"
                    genelements.loc[base, 'eltype'] = EXON
                    genelements.loc[base, 'start'] = base
                    genelements.loc[base, 'finish'] = base + ex_length -1
                    base = base + ex_length
                    logger.info("Gene added")    
                        
                    while (np.random.random_sample() < 0.75):  #25% probability of stopping
                        in_length = round(np.random.normal(450, 100))
                        genelements.loc[base, 'name'] = "intron"
                        genelements.loc[base, 'eltype'] = INTRON
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + in_length -1
                        base = base + in_length;
                        
                        ex_length = round(np.random.lognormal(math.log(250), math.log(1.3))) + 1
                        genelements.loc[base, 'name'] = "exon"
                        genelements.loc[base, 'eltype'] = EXON
                        genelements.loc[base, 'start'] = base
                        genelements.loc[base, 'finish'] = base + ex_length -1
                        base = base + ex_length 
                            
            #final non-coding region
            genelements.loc[base, 'name'] = "noncoding"
            genelements.loc[base, 'eltype'] = NONCOD
            genelements.loc[base, 'start'] = base
            genelements.loc[base, 'finish'] = self.gensize - 1
            logger.info("Chromosome complete!")
            logger.debug(genelements)
            self.genelements = genelements

        elif self.type == "dict":
            pass

        else:
            print("'type' must be 'random' or 'dict'")

    def review(self, item = None):
        """
        allows user to inspect chromosome settings, including:
        mutation rates, genomic element types, and genomic elements
        also creates interactive plot of the chromosome, so the user can inspect it
        """

        if item == "mutations":
            pass
            #print("Mutations:\n" self.mutdict)
        elif item == "eltypes":
            pass
            #print("Genomic Element Types:\n" self.eldict)
        elif item == "elements":
            print("Genomic elements:\n", self.genelements)
        elif item == "chromosome":
            #Make the `rectangles` dataframe for plotting
            eltype = []
            startbase = []
            endbase = []
            y1 = []
            y2 = []
            for index, row in self.genelements.iterrows():
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
            genecount = (rectangles["Element Type"]=='noncoding').sum()-1
            exoncount = (rectangles["Element Type"]=='exon').sum()
            introncount = (rectangles["Element Type"]=='intron').sum()
            totexonlength = 0
            totintronlength = 0
            for index, row in rectangles.iterrows():
                if row["Element Type"]=='exon':
                    totexonlength += (row["x2"]-row["x1"])
                elif row["Element Type"]=='intron':
                    totintronlength += (row["x2"]-row["x1"])


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
            for index, row in self.genelements.iterrows():
                eltype.append(row['name'])
                startbase.append(row['start'])
                endbase.append(row['finish'])
                y1.append(0)
                y2.append(1)

            chromcoords = list(zip(eltype, startbase, endbase, y1, y2))
            rectangles = pd.DataFrame(chromcoords, columns = ['Element Type', 'x1', 'x2', 'y1', 'y2'])

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
                color=alt.Color('Element Type', 
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

            IPython.display.display_html(zoom)

        else:
            logger.info("Please enter a valid item to review. Options include:\n"
                "'mutations'\n'eltypes'\n'elements'\n''chromosome'\n'interractive'")


if __name__ == "__main__":

    # generate random chromosome
    init_chromosome = Chromosome()
    Chromosome.make_rand(init_chromosome)
    Chromosome.review(init_chromosome, item = "elements")
    Chromosome.review(init_chromosome, item = "chromosome")
    Chromosome.review(init_chromosome, item = "interactive")