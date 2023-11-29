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

from typing import Mapping, Sequence, Optional, ClassVar
from dataclasses import dataclass, field
import numpy as np
import scipy.stats as stats
from scipy.stats import rv_continuous
import toyplot
from loguru import logger


DISTOPTS = ['f', 'g', 'e', 'n', 'w', 's']
BAD_DIST_TYPE = """\
Distribution type options: \n"
    "'f' = fixed fitness effect (uniform)\n"
    "'g' = gamma distribution\n"
    "'e' = exponential distribution\n"
    "'n' = normal distribution\n"
    "'w' = Weibull distribution\n"
    "'s' = Script-based distribution"
"""

@dataclass(eq=False)
class MutationType:
    """MutationType Class for fitness effects.

    This is used to model fitness effects of mutations as random
    variables drawn from a probability distribution in SLiM. The
    distributions can be visualized and summarized using convenience
    functions here. The code to define the object in SLIM can be
    generated with the `.to_slim()` function call.

    In contrast to MutationTypes in SLiM, this class includes the
    additional parameter `diff_expr` which applies easily to writing
    alternation of generations (complex life cycle) scripts in shadie.
    With this you can specify whether a mutation affects fitness in
    only the diploid, only the haploid, or both life stages.

    Parameters
    ----------
    dominance: float
        Dominance coefficient scales the effect of the selection
        coefficient on fitness in heterozygotes.
    distribution: str
        A string to represent a distribution of fitness effects.
        Examples are: 'f', 'g', 'n', 'e', 'w', 's'.
    parameters: Sequence[float]
        One or more parameters of the specified distribution.
    affects_diploid: bool
        Affects diploid fitness. Default=True.
    affects_haploid: bool
        Affects haploid fitness. Default=True.
    idx: int
        Set a unique integer index for this MutationType.

    # diff_expr: Optional[str]
    # ...

    Example
    -------
    Create a MutationType object using the `shadie.mtype` function.
    >>> mut0 = shadie.mtype(0.5, 'e', (2.5))
    >>> mut0.print_summary()
    """
    dominance: float
    """: Scales the effect of selection coefficient on diploid heterozygotes."""
    distribution: str
    """: String name for fitness effect distribution: ['f', 'g', 'n', 'w', 'e']"""
    params: Sequence[float]
    """: Parameters of the distribution of fitness effects."""
    affects_diploid: bool
    """: Mutation will affect diploid fitness"""
    affects_haploid: bool
    """: Mutation will affect haploid fitness"""
    convert_to_substitution: bool = False
    """: Determines whether a Mutation is converted to a Substitution object in SLiM"""

    # Class variable
    idx: ClassVar[int] = field(default=0, init=True)

    # Parameters auto-filled and not entered by users
    """: Unique index of this mutation type."""
    _dist: rv_continuous = field(default=None, init=False, repr=False)
    """: Scipy stats distribution of the .distribution name type."""
    _params: Mapping[str, float] = field(default=dict, init=False, repr=False)
    """: A dict mapping parameter names to values for a distribution."""
    _sign: int = field(default=1, init=False, repr=False)
    """: A value of 1 or -1 depending whether dist allows negative."""

    _force_idx: Optional[int] = field(default=None, init=True, repr=False)
    """Private variable to override auto setting of class variable idx."""

    def __post_init__(self):
        """Set the _dist distribution using the user params."""
        # advance counter for of MutationTypes in existence.
        MutationType.idx += 1
        # set to user entered idx else to Mutation counter value
        self.idx = MutationType.idx if self._force_idx is None else self._force_idx

        # user-entered param can be one or multiple values -> Tuple
        self.params = (
            (self.params,) if isinstance(self.params, float)
            else tuple(self.params))

        # parse user entered 'distribution' and 'params' to
        if self.distribution == "f":
            self._dist = stats.uniform
            self._sign = 1
            self._params = {"loc": self.params[0], "scale": 1e-9}
        elif self.distribution == "n":
            self._dist = stats.norm
            self._sign = 1
            self._params = {"loc": self.params[0], "scale": self.params[1]}
        elif self.distribution == "g":
            self._dist = stats.gamma
            self._sign = -1 if self.params[0] < 0 else 1
            self._params = {"a": abs(self.params[0]), "scale": self.params[1]}
        elif self.distribution == "e":
            self._dist = stats.expon
            self._sign = -1 if self.params[0] < 0 else 1
            self._params = {"loc": 0.0, "scale": abs(self.params[0])}
        else:
            raise ValueError(f"distribution {self.distribution} not recognized.")

    def __repr__(self):
        """Return a string reprentation of the object."""
        value = (
            f"MutationType({self.name}, {self.dominance}, "
            f"{self.distribution}, {self.params},"
            f"{self.affects_diploid}, {self.affects_haploid})"
        )
        return value

    ###################################################################
    # Properties called as .name() or .is_coding()
    ###################################################################
    @property
    def name(self) -> str:
        """Returns a unique name for this MutationType.
s
        When called from within a chromosome object MutationTypes
        are renamed with unique numbering.
        """
        return f"m{self.idx}"

    @property
    def is_coding(self) -> bool:
        """Returns True if MutationType has non-zero fitness effects."""
        if (self.distribution == 'f') and (self.params[0] == 0.):
            return False
        return True

    ###################################################################
    # Main function call for writing to SLIM
    ###################################################################
    def to_slim(self, nuc=False) -> str:
        """Returns the SLIM command to Initialize the MutationType.

        This function is called from within a shadie simulation setup
        by accessing this MutationType object from a Chromosome.
        """
        inner = f"'{self.name}', {self.dominance}, '{self.distribution}', "
        inner += ", ".join(map(str, self.params))
        if nuc:
            initialize = f"initializeMutationTypeNuc({inner});"
        else:
            initialize = f"initializeMutationType({inner});"

        if self.convert_to_substitution:
            convert = f"\n{self.name}.convertToSubstitution = T;"
            to_slim_string = initialize + convert
        else:
            to_slim_string = initialize

        return to_slim_string

    ###################################################################
    # Extract stats from distribution.
    ###################################################################
    @property
    def mean(self) -> float:
        """Return the mean selection coefficient from the distribution"""
        return self._dist.mean(**self._params) * self._sign

    @property
    def std(self) -> float:
        """Return the std of the selection coefficient from the distribution"""
        return self._dist.std(**self._params)

    @property
    def min(self) -> float:
        """Return the min value in the 99% CI of the distribution"""
        interval = self._dist.interval(0.995, **self._params)
        if self._sign < 0:
            return interval[1] * self._sign
        return interval[0]

    @property
    def max(self) -> float:
        """Return the max value in the 99% CI of the distribution"""
        interval = self._dist.interval(0.995, **self._params)
        if self._sign < 0:
            return interval[0] * self._sign
        return interval[1]

    def print_summary(self) -> None:
        """Prints a summary statement about the distribution."""
        print(
            '\033[1m' + "Mutation Type" + '\033[0m' + "\n"
            f"idx: {self.name}\n"
            f"dominance coefficient: {self.dominance}\n"
            f"distribution: {self.distribution}\n"
            f"distribution parameters: {self._params}\n"
            f"distribution_mean: {self.mean:.4f}\n"
            f"distribution_std: {self.std:.4f}\n"
            f"affects diploid fitness: {self.affects_diploid}\n"
            f"affects haploid fitness: {self.affects_haploid}\n"
            f"convert to substitution: {self.convert_to_substitution}\n"
            # f"differential_expr: {self._expr}\n"
        )

    def _draw_hists(self, axes):
        """
        ...
        """
        xpoints = np.linspace(self.min, self.max, 100)
        yvalues = self._dist.pdf(xpoints * self._sign, **self._params)
        mark = axes.fill(
            xpoints, yvalues,
            style={
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
        if self.distribution == "f":
            axes.x.domain.max = self.mean + 1
            axes.x.domain.min = self.mean - 1
        return mark

    def draw(self, axes=None, show_mean=True, **kwargs):
        """Return a toyplot histogram of the selection coefficients."""
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
        self.print_summary()
        return self.draw()


def mtype(
    dominance: float,
    distribution: str,
    params: Sequence[float],
    affects_diploid: bool = True,
    affects_haploid: bool = True,
    convert_to_substitution: bool = False,
    force_idx: Optional[int] = None,
) -> MutationType:
    """MutationType class constructor function.

    Returns a MutationType instance for the selected 'distribution'
    type and its parameter settings.

    Parameters
    ----------
    ...
    """
    # group params as a tuple
    if isinstance(params, float):
        params = (params,)

    # raise exception if bad type
    if distribution not in DISTOPTS:
        logger.info(BAD_DIST_TYPE)
        raise ValueError(f"distribution type must be one of {DISTOPTS}")

    # group all args into a dict
    kwargs = dict(
        dominance=dominance,
        distribution=distribution,
        params=params,
        affects_diploid=affects_diploid,
        affects_haploid=affects_haploid,
        convert_to_substitution=convert_to_substitution,
        _force_idx=force_idx,
    )

    if distribution == 'f':
        assert len(params) == 1, "fixed (uniform) dist requires one param"
        return MutationType(**kwargs)

    if distribution == 'n':
        assert len(params) == 2, "normal dist requires two params (loc, scale)"
        return MutationType(**kwargs)

    if distribution == 'g':
        assert len(params) == 2, "gamma dist requires two params (a, scale)"
        return MutationType(**kwargs)

    if distribution == 'e':
        assert len(params) == 1, "exponential dist requires one param"
        return MutationType(**kwargs)
    raise ValueError(f"distribution {distribution} not recognized.")


if __name__ == "__main__":

    m = mtype(0.5, 'n', [0.0, 1.0])
    m.print_summary()

    m = mtype(0.5, 'f', 0.1)
    m.print_summary()

    m = mtype(0.5, 'g', (2.0, 0.1))
    m.print_summary()

    m = mtype(0.5, 'e', (2.5))
    m.print_summary()

    m = mtype(0.5, 'e', (2.5), force_idx=99)
    m.print_summary()

    m = mtype(0.5, 'e', (2.5), force_idx=0)
    m.print_summary()

    m = mtype(0.5, 'g', (2, 0.1), force_idx=0, affects_diploid=False, convert_to_substitution= True)
    m.print_summary()

    print(m.to_slim())

    #print(m.__doc__)

    # mlist = shadie.mlist(
    #     #shadie.mtype(0.5, 'f', 0.0),
    #     # <<<<<<< HEAD
    #     shadie.mtype(0.5, 'f', 0.1, True, True),
    #     shadie.mtype(0.5, 'n', (0.5, 0.25)),
    #     shadie.mtype(0.5, 'g', (2.0, 0.1)),
    #     shadie.mtype(0.1, 'e', 2.5),
    #     # =======
    #     # SLIM4 EDITION
    #     # shadie.mtype(0.5, 'f', 0.1),
    #     # shadie.mtype(0.5, 'n', 0.5, 0.25),
    #     # shadie.mtype(0.5, 'g', 2.0, 0.1),
    #     # shadie.mtype(0.1, 'e', 2.5, diffexpr = "diploid"),
    #     # >>>>>>> 1283be0eadbb91c01407cc5c09497b96f6a762fa
    # )
    # for muta in mlist:
    #     muta.summary()
    # #m2 = shadie.mtype(0.5, 'f', 0)

    # print(mlist)
    # test= []
    # expr = []
    # for mut in mlist:
    #     test.append(mut.name)
    #     expr.append(mut._expr)
    # print(test)
    # print(mlist[0].to_slim())

    # mut2 = shadie.mtype(0.5, 'f', 0.1, diffexpr="diploid")
    # print(mut2._expr)
