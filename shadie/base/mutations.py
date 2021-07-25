#!/usr/bin/env python

"""
Allows user to create MutationType and MutationList instances

SHADIE usage:
-------------
MutationTypes are used in shadie to describe the regions of a chromosome
and the type of fitness effects that mutations to these regions will 
cause.

SHADIE example:
---------------
mlist = shadie.mlist(
    shadie.mtype(0.5, 'f', 0.1),
    shadie.mtype(0.5, 'n', 0.5, 0.25),
    shadie.mtype(0.5, 'g', 2.0, 0.1),
    shadie.mtype(0.1, 'e', 2.5),
)
for muta in mlist:
    muta.summary()
    print(muta.to_slim())

SLIM example:
-------------
initializeMutationType(1, 0.5, "f", 0);
initializeMutationType(2, 0.5, "f", 0);
initializeMutationType(3, 0.1, "g", -0.03, 0.2);
initializeMutationType(4, 0.8, "e", 0.1);
"""


import numpy as np
import scipy.stats as stats
import toyplot
from loguru import logger


DISTOPTS = ['f', 'g', 'e', 'n', 'w', 's']


class MutationTypeBase:
    """
    Generic mutation class instance to be inherited by a specific type.
    """
    idx = 0

    def __init__(
        self,
        dominance: float,
        distribution: str,
        *params: float
        ):

        MutationTypeBase.idx += 1
        self.idx = MutationTypeBase.idx
        self.name = f"m{self.idx}"

        self.dom = dominance
        self.dist = distribution
        self.distparams = params
        if self.dist == 'f':
            if self.distparams == (0,):
                self.coding = 0
            else: 
                self.coding = 1
        else:
            self.coding =  1

        # values overwritten by inherited classes
        self._dist = stats.norm
        self._params = {'loc': 0, 'scale': 1}
        self._neg = 1

    def __repr__(self):
        view = [self.name, self.dom, self.dist, tuple(self.distparams)]
        return f"<MutationType: {', '.join(map(str, view))}>"

    def to_slim(self, nuc=False):
        """
        Returns the SLIM command to Initialize the MutationType
        """
        inner = f"'{self.name}', {self.dom}, '{self.dist}', "
        inner += ", ".join(map(str, self.distparams))
        if nuc:
            return (f"initializeMutationTypeNuc({inner});\n    "
                f"{self.name}.convertToSubstitution = T;")
        return (f"initializeMutationType({inner});\n    "
                f"{self.name}.convertToSubstitution = T;")

    @property
    def mean(self):
        "the mean selection coefficient from the distribution"
        return self._dist.mean(**self._params) * self._neg

    @property
    def std(self):
        "the std of the selection coefficient from the distribution"        
        return self._dist.std(**self._params)

    @property
    def min(self):
        "the min value in the 99% CI of the distribution"                
        interval = self._dist.interval(0.99, **self._params)
        if self._neg < 0:
            return interval[1] * self._neg
        return interval[0]

    @property
    def max(self):
        "the max value in the 99% CI of the distribution"
        interval = self._dist.interval(0.99, **self._params)
        if self._neg < 0:
            return interval[0] * self._neg
        return interval[1]


    def summary(self):
        """
        prints a summary statement about the distribution
        """
        print(
            '\033[1m' + "Mutation Type" + '\033[0m' + "\n"
            f"idx: {self.name}\n"
            f"dominance coefficient: {self.dom}\n"
            f"distribution: {self.dist}\n"
            f"distribution parameters: {self.distparams}\n"
            f"distribution_mean: {self.mean:.4f}\n"
            f"distribution_std: {self.std:.4f}\n"
        )


    def _draw_hists(self, axes):
        """
        ...
        """
        xpoints = np.linspace(self.min, self.max, 100)
        yvalues = self._dist.pdf(xpoints * self._neg, **self._params)
        mark = axes.fill(
            xpoints, yvalues, 
            style = {
                "stroke": toyplot.color.Palette()[0], 
                "stroke-width": 1.5,
                "fill": toyplot.color.Palette()[0],             
                "fill-opacity": 0.33,
            },           
        )
        axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
        axes.x.ticks.show = True
        return mark


    def _draw_means(self, axes):
        """
        ...
        """
        # add mean line to the axes
        mark = axes.vlines(
            self.mean, 
            style={
                "stroke": toyplot.color.Palette()[1], 
                "stroke-width": 1.5,
            }
        )
        if self.dist == "f":
            axes.x.domain.max = self.mean + 1
            axes.x.domain.min = self.mean - 1
        return mark


    def draw(self, axes=None, show_mean=True, **kwargs):
        """
        Returns a toyplot histogram of the selection coefficients.
        """
        # create new axes if none was provided
        canvas = toyplot.Canvas(
            width=kwargs.get('width', 300), 
            height=kwargs.get('height', 225),
        )
        axes = canvas.cartesian(
            xlabel="fitness effect",
            yshow=False,
        )
        marks = []
        marks.append(self._draw_hists(axes))
        if show_mean:
            marks.append(self._draw_means(axes))
        return canvas, axes, marks


    def inspect(self):
        """
        ...
        """
        self.summary()
        return self.draw()


class MutationNormal(MutationTypeBase):
    def __init__(self, dominance, shape, scale):
        super().__init__(dominance, 'n', shape, scale)
        self._dist = stats.norm
        self._neg = 1
        self._params = {"loc": shape, "scale": scale}


class MutationGamma(MutationTypeBase):
    def __init__(self, dominance, shape, scale):
        super().__init__(dominance, 'g', shape, scale)
        self._dist = stats.gamma
        self._neg = (-1 if shape < 0 else 1)
        self._params = {"a": abs(shape), "scale": scale}


class MutationFixed(MutationTypeBase):
    def __init__(self, dominance, shape):
        super().__init__(dominance, 'f', shape)
        self._dist = stats.uniform
        self._neg = 1
        self._params = {"loc": shape, "scale": 1e-9}


class MutationExponential(MutationTypeBase):
    def __init__(self, dominance, scale):
        super().__init__(dominance, 'e', scale)
        self._dist = stats.expon
        self._neg = (-1 if scale < 0 else 1)
        self._params = {"loc": 0.0, "scale": abs(scale)}


# class MutationWeibull(MutationTypeBase):
#     def __init__(self, dominance, shape):
#         super().__init__(dominance, 'e', scale)
#         self._dist = stats.expon
#         self._neg = (-1 if scale < 0 else 1)
#         self._params = {"loc": 0.0, "scale": abs(scale)}



def mtype(dominance, distribution, *params):
    """
    MutationType class constructor. Returns a MutationType subclass 
    instance for the selected 'distribution' type.

    Parameters:
    -----------
    ...
    """
    if distribution not in DISTOPTS:
        logger.info("Distribution type options: \n"
            "'f' = fixed fitness effect\n"
            "'g' = gamma distribution\n"
            "'e' = exponential distribution\n"
            "'n' = normal distirbution\n"
            "'w' = Weibull distribution\n"
            "'s' = Script-based distribution")
        raise ValueError(f"distribution type must be one of {DISTOPTS}")

    if distribution == 'n':
        assert len(params) == 2, "normal dist requires two params (loc, scale)"
        return MutationNormal(dominance, *params)
    if distribution == 'g':
        assert len(params) == 2, "gamma dist requires two params (a, scale)"
        return MutationGamma(dominance, *params)
    if distribution == 'f':
        assert len(params) == 1, "fixed (uniform) dist requires one param"
        return MutationFixed(dominance, *params)
    if distribution == 'e':
        assert len(params) == 1, "exponential dist requires one param"
        return MutationExponential(dominance, *params)
    raise NotImplementedError("distribution not yet supported; TODO")



class MutationList(list):
    """
    Creates a list of mutations and the mutationdict object for Shadie
    """
    def __init__(self, *mutationtypes):
        super().__init__(mutationtypes)

        self.names = [i.name for i in self]
        self.dict = {i.name: i for i in self}
        self.slim_dict = {}        

        # build a string representation of muttype for slim
        for mut in self:
            params = ", ".join(map(str, mut.distparams))
            srep = f"'{mut.name}', {mut.dom}, '{mut.dist}', {params}"
            self.slim_dict[mut.name] = srep

    def __repr__(self):
        return f"<MutationList: {self.names}>"

    @property
    def min(self):
        "min value in the 99% CI of all MutationType distributions"
        return min([i.min for i in self])

    @property
    def max(self):
        "max value in the 99% CI of all MutationType distributions"
        return max([i.max for i in self])


    def draw(self, **kwargs):
        """
        Returns a toyplot histogram of the selection coefficients.
        """
        canvas, axes, marks = self[0].draw(show_mean=False, **kwargs)
        for mut in self[1:]:
            marks.append(mut._draw_hists(axes=axes))
        return canvas, axes, marks


if __name__ == "__main__":

    # generate random chromosome
    import shadie

    mlist = shadie.mlist(
        shadie.mtype(0.5, 'f', 0.0),
        shadie.mtype(0.5, 'f', 0.1),
        shadie.mtype(0.5, 'n', 0.5, 0.25),
        shadie.mtype(0.5, 'g', 2.0, 0.1),
        shadie.mtype(0.1, 'e', 2.5),
    )
    for muta in mlist:
        muta.summary()
    m2 = shadie.mtype(0.5, 'f', 0)

    print(mlist)
    test= []
    for mut in mlist:
        test.append(mut.name)
    print(test)
    print(mlist[0].to_slim())
