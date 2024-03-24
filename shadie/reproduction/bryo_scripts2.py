#!/usr/bin/env python

"""
...
"""


FUNCTION_MAKE_HOMOSPORES = """
function (void)make_homospores(object<Individual>$ ind, integer$ reps) {
    for (rep in 1:reps){
        breaks1 = sim.chromosome.drawBreakpoints(ind);
        breaks2 = sim.chromosome.drawBreakpoints(ind);
        p0.addRecombinant(ind.genome1, ind.genome2, breaks1, NULL, NULL, NULL).tag=3;
        p0.addRecombinant(ind.genome2, ind.genome1, breaks1, NULL, NULL, NULL).tag=3;
        p0.addRecombinant(ind.genome1, ind.genome2, breaks2, NULL, NULL, NULL).tag=3;
        p0.addRecombinant(ind.genome2, ind.genome1, breaks2, NULL, NULL, NULL).tag=3;
    }
}
"""