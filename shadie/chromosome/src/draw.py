#!/usr/bin/env python

"""
Inspect tools for visualizing chromosome structure with alt.
"""

from typing import Optional
import pandas as pd
import altair as alt
# from shadie.chromosome.src.base_class import ChromosomeBase


def draw_altair_chrom_canvas(chrom: 'ChromosomeBase', width: int=700):
    """Return an altair drawing of the chromosome."""

    # collect data from chromosome
    data = chrom.data.copy()
    data['category'] = chrom.data.name
    data['altname'] = chrom.data.script.apply(lambda x: x.altname)
    data['length'] = chrom.data.end - chrom.data.start

    # if name was empty then infer type from coding status
    for idx in data.index:
        if pd.isna(data.category[idx]):
            if data.coding[idx]:
                data.loc[idx, "category"] = "exon"
            elif ("intron" in data.name[idx]) or ("intron" in data.script.altname):
                data.loc[idx, "category"] = "intron"
            else:
                data.loc[idx, "category"] = "noncds"

    # subset to the data used for plotting
    genome = data[["eltype", "category", "altname", "length"]].copy()
    genome["x1"] = data.start
    genome["x2"] = data.end
    genome["y1"] = 0
    genome["y2"] = 1

    # set the colors
    cmap_in = ['mediumaquamarine', 'palegreen','olivedrab', 'darkgreen', 'limegreen']
    cmap_ex = ['cornflowerblue', 'mediumblue','dodgerblue', 'darkslateblue', 'skyblue']
    cmap_nc = ['lemonchiffon', 'gold', 'orange',  'yellow', 'khaki']

    # count regions for applying colors
    n_ex = genome.loc[genome.category=='exon'].altname.nunique()
    t_ex = genome.loc[genome.category=='exon'].altname.unique().tolist()
    n_in = genome.loc[genome.category=='intron'].altname.nunique()
    t_in = genome.loc[genome.category=='intron'].altname.unique().tolist()
    n_nc = genome.loc[genome.category=='noncds'].altname.nunique()
    t_nc = genome.loc[genome.category=='noncds'].altname.unique().tolist()

    # colors in regions
    dom = t_ex + t_in + t_nc
    rng = cmap_ex[:n_ex] + cmap_in[:n_in] + cmap_nc[:n_nc]

    # create a chromosome composed of rectangles
    ichrom = (alt.Chart(genome)
        .mark_rect()
            .encode(
                x=alt.X('x1:Q', axis=alt.Axis(title='Base Pairs')),
                y=alt.Y('y1:Q', axis=None),
                x2='x2:Q',
                y2='y2:Q', 
                color=alt.Color(
                    'altname:N', 
                    scale=alt.Scale(domain=dom, range=rng)
                ),
                tooltip=[
                    alt.Tooltip('eltype', title='Element Type'),
                    alt.Tooltip('altname', title='Name'),
                    alt.Tooltip('length', title='Length'),
                    alt.Tooltip('x1', title='Start'),
                    alt.Tooltip('x2', title='Stop'),
                ]
            ).properties(width=width)
        )
    return ichrom


def draw_altair_chrom_canvas_interactive(
    chrom: 'ChromosomeBase',
    width: int=700, 
    outfile:Optional[str]=None,
    ):
    """Return an interactive altair visualization of the chromosome.
    """
    # create a brush to select interval
    mark = alt.BrushConfig(fill='red', fillOpacity=0.700, stroke="black")
    brush = alt.selection_interval(encodings=['x'], mark=mark)
    
    # get the canvas drawing    
    ichrom = draw_altair_chrom_canvas(chrom, width=width)

    # create an interactive composite plot with two views of ichrom
    # the first view is 2X as tall and is scaled by the brush selector
    view1 = ichrom.encode(
        alt.X('x1:Q', title=None, scale=alt.Scale(domain=brush))
        ).properties(height=80)
    # the second view has the brush selector active.
    view2 = ichrom.add_selection(brush).properties(height=40)
    zoom = alt.vconcat(view1, view2, data=ichrom.data)

    # optionally write to disk
    if outfile:
        zoom.save(outfile.strip('.html') + '.html')
    return zoom


def toyplot(self):
    "Makes static toyplot"
    self.inspect()
    eltype = []
    startbase = []
    endbase = []
    y1 = []
    y2 = []
    for index, row in self.genome.iterrows():
        eltype.append(row['category'])
        startbase.append(row['x1'])
        endbase.append(row['x2'])
        y1.append(0)
        y2.append(1)

    chromcoords = list(zip(eltype, startbase, endbase, y1, y2))
    rectangles = pd.DataFrame(chromcoords, columns = [
        'Element Type', 'x1', 'x2', 'y1', 'y2'])
    self.rectangles = rectangles

    #define colors for each element
    color = []
    for index, row in rectangles.iterrows():
        if row["Element Type"] == "noncds":
            color.append("lemonchiffon")
        elif row["Element Type"] == "exon":
            color.append("royalblue")
        elif row["Element Type"] == "intron":
            color.append("mediumaquamarine")

    #assign colors to rectangles dataframe
    rectangles.insert(5, "color", color)

    #Calculate Stats
    self.genecount = max(1, (rectangles["Element Type"]=='noncds').sum()-1)
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


def review(self, item=None):
    """Return a summary and visualization of the chrosome structure.
    
    Select an item to review among the following options:
        mutations: ...
        eltypes: ...
        elements: ...
        chromosome: ...
        interractive: ...
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
        self.alt()
        print("Interactive alt chromosome map:")
        IPython.display.display_html(self.zoom)
        
    else:
        logger.info("Please enter a valid item to review. Options include:\n"
            "'mutations'\n'eltypes'\n'elements'\n''chromosome'\n'interractive'")


if __name__ == "__main__":

    import shadie
    test_chrom = shadie.chromosome.default()
    print(draw_altair_chrom_canvas_interactive(test_chrom))
