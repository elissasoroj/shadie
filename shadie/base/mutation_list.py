#!/usr/bin/env python

"""MutationList objects are used to...

"""

from typing import List, Dict
from shadie.base.mutations import MutationType


class MutationList(list):
    """A List child class to hold multiple MutationType objects.

    This object has additional attributes and methods for summarizing
    or visualizing the effects of multiple MutationTypes, and for
    writing them to SLiM code format.
    """
    def __init__(self, *mutationtypes):
        super().__init__(mutationtypes)

    def __repr__(self):
        return f"<MutationList: {self.names}>"

    @property
    def names(self) -> List[str]:
        """Return a list of the names of Mutations in list."""
        return [i.name for i in self]

    @property
    def dict(self) -> Dict[str, MutationType]:
        """Return a dict mapping MutationType names to MutationType objects."""
        return {i.name: i for i in self}

    @property
    def slim_dict(self) -> Dict[str, str]:
        """Return dict mapping MutationType names to SLiM code for MutationType."""
        slim_dict = {}
        for mut in self:
            # build a string representation of muttype for slim
            params = ", ".join(map(str, mut.params))
            srep = f"'{mut.name}', {mut.dominance}, '{mut.distribution}', {params}"
            slim_dict[mut.name] = srep
        return slim_dict

    @property
    def min(self):
        """min value in the 99% CI of all MutationType distributions"""
        return min([i.min for i in self])

    @property
    def max(self):
        """max value in the 99% CI of all MutationType distributions"""
        return max([i.max for i in self])

    def draw(self, **kwargs):
        """Returns a toyplot histogram of the selection coefficients."""
        canvas, axes, marks = self[0].draw(show_mean=False, **kwargs)
        for mut in self[1:]:
            marks.append(mut._draw_hists(axes=axes))
        return canvas, axes, marks


if __name__ == "__main__":
    pass
