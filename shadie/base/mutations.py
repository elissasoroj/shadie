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


from typing import Iterable, Union, List, Dict
from dataclasses import dataclass, field
import numpy as np
import scipy.stats as stats
import toyplot
from loguru import logger


DISTOPTS = ['f', 'g', 'e', 'n', 'w', 's']


@dataclass
class MutationTypeBase:
    """Mutation Base Class to be inherited by more specific types.

    The functions of this object require additional parameters that
    are only created in the superclasses.

    Parameters
    ----------
    dominance: float
        Dominance coefficient scales the effect of the selection 
        coefficient on fitness in heterozygotes.
    sporophyte_phenotype: bool
        If True then the selection coefficient affects fitness during
        the sporophyte generation.
    gametophyte_phenotype: bool
        If True then the selection coefficient affects fitness during
        the gametophyte generation.
    """
    idx: int = field(default=0, init=False)
    dominance: float
    sporophyte_phenotype: bool
    gametophyte_phenotype: bool

    distribution: str = field(default='f', init=False)
    params: List[float] = field(default=0.0, init=False)

    def __repr__(self):
        """Return a string reprentation of the object."""
        value = (
            f"MutationType({self.name}, {self.dominance}, "
            f"{self.distribution}, {self.params})"
        )
        return value

    @property
    def name(self):
        """Returns a unique name for this MutationType.

        When called from within a chromosome object MutationTypes
        are renamed with unique numbering.
        """
        return f"m{self.idx}"

    @property
    def coding(self):
        """Returns True if MutationType has non-zero fitness effects."""
        if (self.distribution == 'f') and (self.params[0] == 0.):
            return False
        return True

    def to_slim(self, nuc=False):
        """Returns the SLIM command to Initialize the MutationType.
        
        This function is called from within a shadie simulation setup
        by accessing this MutationType object from a Chromosome.
        """
        inner = f"'{self.name}', {self.dominance}, '{self.distribution}', "
        inner += ", ".join(map(str, self.params))
        if nuc:
            return f"initializeMutationTypeNuc({inner});"
        return f"initializeMutationType({inner});"


@dataclass
class MutationTypeParameterized(MutationTypeBase):
    """Parent class of MutationType classes with stat distributions.

    Parameters
    ----------
    distribution: str
        The statistical distribution to draw selection coefficients
        from. Supported options are f,g,e,n,w,s.
    params
        Parameters of the specified distribution.
    """
    distribution: str
    params: List[float]
    dist: 'scipy.stats.continuous_distns' = field(init=False)
    params_dict: Dict[str, float] = field(init=False)
    sign: int = field(init=False)

    def __post_init__(self):
        if isinstance(self.params, float):
            self.params = (self.params,)
        if self.distribution == "f":
            self.dist = stats.uniform
            self.sign = 1
            self.params_dict = {"loc": self.params[0], "scale": 1e-9}
        elif self.distribution == "n":
            self.dist = stats.norm
            self.sign = 1
            self.params_dict = {"loc": self.params[0], "scale": self.params[1]}
        elif self.distribution == "g":
            self.dist = stats.gamma
            self.sign = (-1 if self.params[0] < 0 else 1)
            self.params_dict = {"a": abs(self.params[0]), "scale": self.params[1]}
        elif self.distribution == "e":
            self.dist = stats.expon
            self.sign = (-1 if self.params[0] < 0 else 1)
            self.params_dict = {"loc": 0.0, "scale": abs(self.params[0])}
        else:
            raise ValueError(f"distribution {self.distribution} not recognized.")

    @property
    def mean(self):
        """Return the mean selection coefficient from the distribution"""
        return self.dist.mean(**self.params_dict) * self.sign

    @property
    def std(self):
        """Return the std of the selection coefficient from the distribution"""
        return self.dist.std(**self.params_dict)

    @property
    def min(self):
        """Return the min value in the 99% CI of the distribution"""
        interval = self.dist.interval(0.99, **self.params_dict)
        if self.sign < 0:
            return interval[1] * self.sign
        return interval[0]

    @property
    def max(self):
        """Return the max value in the 99% CI of the distribution"""
        interval = self.dist.interval(0.99, **self.params_dict)
        if self.sign < 0:
            return interval[0] * self.sign
        return interval[1]

    def summary(self):
        """Prints a summary statement about the distribution."""
        print(
            '\033[1m' + "Mutation Type" + '\033[0m' + "\n"
            f"idx: {self.name}\n"
            f"dominance coefficient: {self.dominance}\n"
            f"distribution: {self.distribution}\n"
            f"distribution parameters: {self.params}\n"
            f"distribution_mean: {self.mean:.4f}\n"
            f"distribution_std: {self.std:.4f}\n"
        )

    def _draw_hists(self, axes):
        """
        ...
        """
        xpoints = np.linspace(self.min, self.max, 100)
        yvalues = self.dist.pdf(xpoints * self.sign, **self.params_dict)
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


def mtype(
    dominance: float, 
    distribution: str, 
    params: Union[float, Iterable[float]],
    sporophyte_phenotype: bool=True, 
    gametophyte_phenotype: bool=True,
    ):
    """
    MutationType class constructor. Returns a MutationType subclass 
    instance for the selected 'distribution' type.

    Parameters:
    -----------
    ...
    """
    if isinstance(params, float):
        params = (params,)

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
        return MutationTypeParameterized(
            dominance, sporophyte_phenotype, gametophyte_phenotype, 'n', params)
    if distribution == 'g':
        assert len(params) == 2, "gamma dist requires two params (a, scale)"
        return MutationTypeParameterized(
            dominance, sporophyte_phenotype, gametophyte_phenotype, 'g', params)
    if distribution == 'f':
        assert len(params) == 1, "fixed (uniform) dist requires one param"
        return MutationTypeParameterized(
            dominance, sporophyte_phenotype, gametophyte_phenotype, 'f', params)
    if distribution == 'e':
        assert len(params) == 1, "exponential dist requires one param"
        return MutationTypeParameterized(
            dominance, sporophyte_phenotype, gametophyte_phenotype, 'e', params)
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
            params = ", ".join(map(str, mut.params))
            srep = f"'{mut.name}', {mut.dominance}, '{mut.distribution}', {params}"
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
        """Returns a toyplot histogram of the selection coefficients.
        """
        canvas, axes, marks = self[0].draw(show_mean=False, **kwargs)
        for mut in self[1:]:
            marks.append(mut._draw_hists(axes=axes))
        return canvas, axes, marks


if __name__ == "__main__":

    # generate random chromosome
    import shadie

    mlist = shadie.mlist(
        #shadie.mtype(0.5, 'f', 0.0),
        shadie.mtype(0.5, 'f', 0.1, True, True),
        shadie.mtype(0.5, 'n', (0.5, 0.25)),
        shadie.mtype(0.5, 'g', (2.0, 0.1)),
        shadie.mtype(0.1, 'e', 2.5),
    )
    for muta in mlist:
        muta.summary()
    #m2 = shadie.mtype(0.5, 'f', 0)

    print(mlist)
    test= []
    for mut in mlist:
        test.append(mut.name)
    print(test)
    print(mlist[0].to_slim())
