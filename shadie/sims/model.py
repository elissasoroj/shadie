#!/usr/bin/env python

"""
A context manager for wrapping the context of a simulation setup.

initialize() {
    initializeSLiMModelType("nonWF");
    defineConstant("K", 500);
    
    initializeMutationType("m1", 0.5, "f", 0.0);
    m1.convertToSubstitution = T;                      // default is T
    
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, 99999);
    initializeMutationRate(1e-7);
    initializeRecombinationRate(1e-8);
}
reproduction() {
    subpop.addCrossed(individual, subpop.sampleIndividuals(1));
}
1 early() {
    sim.addSubpop("p1", 10);
}
early() {
    p1.fitnessScaling = K / p1.individualCount;
}
late() {
    inds = p1.individuals;
    catn(sim.generation + ": " + size(inds) + " (" + max(inds.age) + ")");
}
2000 late() {
    sim.outputFull(ages=T);
}

"""

import subprocess
from typing import Union
from contextlib import AbstractContextManager
from loguru import logger
from shadie.base.mutations import MutationTypeBase
from shadie.base.elements import ElementType
# from shadie.reproduction.reproduction import Reproduction
#from shadie.reproduction import Reproduction

# cannot do both mutationRate and nucleotidebased 

INIT = """
initialize() {{
   
  // model type
  initializeSLiMModelType("nonWF");

  // config
  initializeRecombinationRate({recombination_rate});
  initializeMutationRate({mutation_rate});
  initializeTreeSeq();

  // MutationType init
  {mutations}

  // ElementType init
  {elements}

  // Chromosome (GenomicElement init)
  {chromosome}

  // constants (Ne)
  {constants}

  // extra scripts 
  {scripts}
}}
"""
# --------------------------------------------

REPRO = """
reproduction({population}) {{ //generates offspring
    {scripts}
}}
"""

FIT = """
{idx} fitness({mutation}) //adjusts fitness calculation
{{
    {scripts}

}}
"""

SURV = """
{idx} survival({population}) //implements survivavl adjustments
{{
    {scripts}

}}
"""

# --------------------------------------------

EARLY = """
{time} early() //executes after offspring are generated
{{
  {scripts}
}}
"""


LATE = """
{time} late() //execuutes after selection occurs
{{
  {scripts}
}}
"""


class Model(AbstractContextManager):
    """
    ...

    Parameters
    -----------
    ...

    """
    def __init__(self, ):
        
        # hold script components as a dict until __exit__ converts to a str
        self.script = {}
        self.stdout = ""

        # store chromosome (mut, ele), constants, and populations
        self.chromosome = None
        self.constants = {}
        self.populations = {}
        self.length = {} #length of simulation in generations

        # self.reproduction = Reproduction()


    def __repr__(self):
        return "<shadie.Model ... >"


    def __enter__(self):
        """
        On entry the Class counters of mutation and element types
        is reset to zero.
        """
        MutationTypeBase.idx = 0
        ElementType.idx = 0
        return self


    def __exit__(self, exc_type, exc_value, exc_traceback):
        """
        Fills the script with context-defined content
        """
        # order script keys
        nulls = [i for i in self.script if i[1] is None]
        ikeys = sorted([i for i in self.script if isinstance(i[1], int)])

        # which type is not null or integer?
        lkeys = [i for i in self.script if i not in nulls or ikeys]

        # remove init from the nulls b/c it must come first
        sorted_keys = [nulls.pop(nulls.index(("initialize", None)))]

        # then order keys by: sorted-integers, nulls, others
        sorted_keys += ikeys
        sorted_keys += nulls
        sorted_keys += lkeys

        # compress script to string
        self.script = "\n".join([self.script[i] for i in sorted_keys])

        # attempt to validate script
        self._check_script()
        logger.debug("exiting Model")


    def initialize(
        self,
        chromosome,
        length:int=1000, #length of sim in # of generations
        mut:float=1e-8, 
        recomb:float=1e-9, 
        constants:Union[None, dict]=None,
        scripts:Union[None, list]=None,
        fileout:str="shadie.trees",
        ):
        """
        Initialize a simulation. This fills the SLIM intialize() code
        block with MutationType, ElementType, Element, and other init
        type code.
        """
        logger.debug("initializing Model")
        constants = {} if constants is None else constants
        scripts = [] if scripts is None else scripts
        
        self.script[('initialize', None)] = (
            INIT.format(**{
                "mutation_rate": mut,
                "recombination_rate": recomb,
                "genome_size": chromosome.genome_size,
                "mutations": chromosome.to_slim_mutation_types(),
                "elements": chromosome.to_slim_element_types(),
                "chromosome": chromosome.to_slim_elements(),
                "constants": "\n  ".join([
                    f"defineConstant('{key}', {val});" for key, val
                    in constants.items()
                ]),
                "scripts": "\n  ".join([i.strip(";") + ";" for i in scripts]),
            })
        )


    def repro(self, population:Union[str, None], scripts:Union[str, list]
        ):
        """
        Add reproduction block code here.
        """
        logger.debug("Reproduction Block")

        # compress list of scripts into a string
        if isinstance(scripts, list):
            scripts = "\n  ".join([i.strip(';') + ';' for i in scripts])

        # population as str or empty
        pop_str = str(population) if population else ""

        self.script[("reproduction", population)] = (
            REPRO.format(**{'population': pop_str, 'scripts': scripts})
        ).lstrip()


    def early(self, time:Union[int, None], scripts:Union[str, list]):
        """
        Add an event that happens before selection in every generation (early).
        """
        # todo: validate script
        # compress list of scripts into a string
        if isinstance(scripts, list):
            scripts = "\n  ".join([i.strip(';') + ';' for i in scripts])

        # time as int or empty
        time_str = str(time) if time else ""

        # expand EARLY script block
        self.script[("early", time)] = (
            EARLY.format(**{'time': time_str, 'scripts': scripts})
        ).lstrip()


    def fitness(self, mutation:Union[str, None], scripts:Union[str, list], idx:Union[str, None]):
        """
        Add an event that adjusts fitness values before fitness calculation.
        """
        # todo: validate script
        # compress list of scripts into a string
        if isinstance(scripts, list):
            scripts = "\n  ".join([i.strip(';') + ';' for i in scripts])

        # idx as str or empty
        idx_str = str(idx) if idx else ""
        # mutation as str or empty
        mutation_str = str(mutation) if mutation else ""

        # expand FITNESS script block
        self.script[("fitness", mutation)] = (
            FIT.format(**{'idx': idx_str, 'mutation': mutation_str, 'scripts': scripts})
        ).lstrip()


    def survival(self, population:Union[str, None], scripts:Union[str, list], idx:Union[str, None]):
        """
        Add an event that adjusts fitness values before fitness calculation.
        """
        # todo: validate script
        # compress list of scripts into a string
        if isinstance(scripts, list):
            scripts = "\n  ".join([i.strip(';') + ';' for i in scripts])

        # idx as str or empty
        idx_str = str(idx) if idx else ""
        # mutation as str or empty
        population_str = str(population) if population else ""

        # expand FITNESS script block
        self.script[("survival", population)] = (
            SURV.format(**{'idx': idx_str, 'population': population_str, 'scripts': scripts})
        ).lstrip()


    def late(self, time:Union[int, None], scripts:Union[str, list]):
        """
        Add an event that happens after every generation (late) if 
        time is None, or only after a particular generation if time
        is an integer.
        """
        # todo: validate script
        if isinstance(scripts, list):
            scripts = "\n  ".join([i.strip(';') + ';' for i in scripts])

        # time as int or empty
        time_str = str(time) if time else ""

        # expand LATE script block
        self.script[("late", time)] = (
            LATE.format(**{'time': time_str, 'scripts': scripts})
        ).lstrip()


    def custom(self, scripts:str):
        """
        Add custom scripts outside without formatting by shadie. 
        Scripts must be Eidos-formatted 
        """
        self.script[("custom", None)] = (scripts).lstrip()


    def _check_script(self):
        """
        TODO.
        """
        assert "initialize" in self.script, (
            "You must call initialize() from within Model context")
        # assert "late" in self.script, (
            # "You must call late() from within Model context")


    def write(self, outname):
        """
        Write SLIM script to the outname filepath.
        """
        with open(outname, 'w') as out:
            out.write(self.script)


    def run(self):
        """
        Run SLIM script in subprocess and store STDOUT to self.stdout.
        """
        # write tmp file and call it
        self.write("/tmp/slim.slim")
        cmd = ["slim", "/tmp/slim.slim"]

        # capture stdout
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, _ = proc.communicate()

        # check for errors
        if proc.returncode:
            logger.error(out.decode())
            raise SyntaxError("SLiM3 error")

        # todo: parse stdout to store, and warnings to logger
        self.stdout = out.decode()
        # _, warnings, stdout = out.decode().split("#")
        # self.stdout = stdout
        # logger.warning(warnings.split("//")[0].strip())
        # ...


if __name__ == "__main__":

    import shadie

    with shadie.Model() as model:
        
        # define mutation types
        m0 = shadie.mtype(0.5, 'n', 2.0, 1.0)
        m1 = shadie.mtype(0.5, 'g', 3.0, 1.0)
        
        # define elements types
        e0 = shadie.etype([m0, m1], [1, 2])
        e1 = shadie.etype([m1], [1])
        
        # design chromosome of elements
        chrom = shadie.chromosome.random(
            genome_size=20000,
            noncds=e0,
            intron=e0,
            exon=e1,
        )

        print(chrom.data.head())
        print(chrom.mutations)

        # init the model
        model.initialize(chromosome=chrom)
        
        # add reproduction 
        # model.reproduction()
        # model.early(1000, "sim.addSubpop('p1', 1000); //diploid sporophytes")
        # model.fitness("m4", "return 1 + mut.selectionCoeff; //gametophytes have no dominance effects", "s1" )
        # model.custom("s2 fitness(m5) { return 1 + mut.selectionCoeff; //gametophytes have no dominance effects }")

    print(model.script)
