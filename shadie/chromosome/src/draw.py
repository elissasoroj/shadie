#!/usr/bin/env python

"""
Inspect tools for visualizing chromosome structure with alt.
"""

from typing import Optional
import pandas as pd
import altair as alt
import toyplot


def draw_altair_chrom_canvas(chrom: 'ChromosomeBase', width: int=700):
    """Return an altair drawing of the chromosome."""

    # collect data from chromosome
    data = chrom.data.copy()
    data['altname'] = chrom.data.script.apply(lambda x: x.altname)
    data['length'] = chrom.data.end - chrom.data.start + 1

    # if name was empty then infer type from coding status
    for idx in data.index:
        if ("int" in data.name[idx]):
            data.loc[idx, "category"] = "intron"
        elif data.coding[idx] or ("ex" in data.name[idx]):
            data.loc[idx, "category"] = "exon"
        else:
            data.loc[idx, "category"] = "noncds"

    # subset to the data used for plotting
    genome = data[["eltype", "category", "altname", "length"]].copy()
    genome["x1"] = data.start
    genome["x2"] = data.end
    genome["y1"] = 0
    genome["y2"] = 1

    # set the colors
    cmap_in = ['mediumaquamarine','olivedrab', 'limegreen', 'darkgreen','palegreen',]
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


def draw_toyplot_chrom(
    chrom: 'ChromosomeBase', 
    width: int=700,
    axes: Optional['toyplot.coordinates.Cartesian']=None,
    ):
    """Return a toyplot drawing of the chromosome.

    """
    if axes is None:
        canvas = toyplot.Canvas(width, height=150)
        axes = canvas.cartesian()
    else:
        canvas = None

    # get colors for elements types
    elements = pd.factorize(chrom.data.eltype)[0]
    colormap = toyplot.color.brewer.palette("Spectral", count=max(4, len(set(elements))))
    colors = [colormap[i] for i in elements]

    marks = []
    # plot bars for element types
    for idx, pos in enumerate(chrom.data.index):
        dat = chrom.data.loc[pos]
        mark = axes.fill(
            [dat.start, dat.end],
            [0, 0], 
            [1, 1],
            color=colors[idx],
            opacity=0.75,
            style={"stroke": "black", "stroke-opacity": 1.0},
            title=(
                f"name: {chrom.data.loc[pos, 'name']}\n"
                f"interval: ({dat.start}, {dat.end})\n"
                f"ElementType: {dat.eltype}\n" 
                f"coding: {bool(dat.coding)}"
            )
        )

    axes.y.show = False
    axes.x.ticks.show = True
    axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True, count=8)
    return canvas, axes, mark


if __name__ == "__main__":

    import shadie
    c, a, m = shadie.chromosome.default().draw()
    c, a, m = shadie.chromosome.random().draw()    
    iviz = shadie.chromosome.default().inspect()
