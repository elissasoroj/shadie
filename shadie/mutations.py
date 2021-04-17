#!/usr/bin/env python

"""
Base classes for chromosome building
"""

class SubstitutionBase:
    idx = 99

    def __init__(self, dominance, *args):
        SubstitutionBase.idx += 1
        self.idx = type(self).idx
        self.dominance = dominance
        self.tag = "f"
        self.args = args

    def __del__(self):
        type(self).idx -= 1

    def __repr__(self):
        tup = (f"m{self.idx}", self.dominance, self.tag)
        tup += self.args
        return str(tup)


class Neutral(SubstitutionBase):
    def __init__(self, dominance, *args):
        super().__init__(dominance, *args)
        self.tag = "f"


class Synonymous(SubstitutionBase):
    def __init__(self, dominance, *args):
        super().__init__(dominance, *args)
        self.tag = "f"


class Deleterious(SubstitutionBase):
    def __init__(self, dominance, *args):
        super().__init__(dominance, *args)
        self.tag = "g"

    def __repr__(self):
        tup = (f"m{self.idx}", self.dominance, self.tag)
        tup += self.args
        return str(tup)


class Beneficial(SubstitutionBase):
    def __init__(self, dominance, *args):
        super().__init__(dominance, *args)
        self.tag = "e"

    def __repr__(self):
        tup = (f"m{self.idx}", self.dominance, self.tag)
        tup += self.args
        return str(tup)



if __name__ == "__main__":

    print([Neutral(0.5, 0.0) for i in range(5)])
    print([Synonymous(0.5, 0.0) for i in range(5)])
    print([Deleterious(0.5, 0.0, -0.3, 0.2) for i in range(5)])
    print([Beneficial(0.5, 0.1) for i in range(5)])            
