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
from shadie.reproduction.api import ReproductionApi
from shadie.sims.format import (
    format_event_dicts_to_strings,
    INITIALIZE,
    EARLY,
    LATE,
    FITNESS,
    REPRODUCTION,
    SURVIVAL,
    CUSTOM,
)

# cannot do both mutationRate and nucleotidebased 


class Model(AbstractContextManager):
    """
    The core shadie object class for creating SLiM scripts and 
    validating input arguments using a context manager. 
    """
    def __init__(self):       
        # .map holds script components as a dict until converted to
        # a string upon __exit__(). The dict maps events to lists of
        # dicts with information for the event code block.
        # {
        #   'initialize': {'constants': ..., 'mutations': ...}, 
        #   'early': {'time': 10, script: '...', comment: '...'},
        #   'late': {'time': None, script: '...', comment: '...'},
        # }
        self.map = {
            'initialize': [],
            'reproduction': [],
            'early': [],
            'late': [],
            'fitness': [],
            'survival': [],
            'reproduction': [],
            'custom': [],
        }
        self.script = ""
        self.stdout = ""

        # store chromosome (mut, ele), constants, and populations
        self.chromosome = None
        self.constants = {}
        self.populations = {}
        self.length = {} #length of simulation in generations

        self.reproduction = ReproductionApi(self)


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
        Fills the script with context-defined content in a set 
        order and runs checks on the script.
        """

        sorted_keys = [
            'initialize', 'timed', 'reproduction', 'early',
            'fitness', 'survival', 'late', 'custom',
        ]

        # copy map and split timed events to a new key list
        mapped = self.map.copy()
        mapped['timed'] = []
        mapped_keys = list(mapped.keys())
        for key in mapped_keys:
            if 'time' in mapped[key]:
                mapped['timed'].append(mapped[key].pop(key))

        # visit events by ordered key type
        script_chunks = []
        for key in sorted_keys:

            # sort events within key type
            if key == "timed":
                events = sorted(mapped[key], key=lambda x: x['time'])
            elif key == "fitness":
                events = sorted(mapped[key], key=lambda x: x['idx'])
            elif key == "survival":
                events = sorted(mapped[key], key=lambda x: x['idx'])
            else:
                events = mapped[key]

            # string format each event and add to script chunks list
            for event in events:
                event = format_event_dicts_to_strings(event)

                # string formatting for code blocks
                if key == "initialize":
                    script = INITIALIZE.format(**event)

                elif key == "early":
                    script = EARLY.format(**event)
                
                elif key == "late":
                    script = LATE.format(**event)

                elif key == 'fitness':
                    script = FITNESS.format(**event)

                elif key == 'survival':
                    script = SURVIVAL.format(**event)

                elif key == 'custom':
                    script = CUSTOM.format(**event)

                elif key == 'reproduction':
                    script = REPRODUCTION.format(**event)

                else:
                    raise NotImplementedError(f"{key} not supported")
                script_chunks.append(script)

        # collapse into the final script string
        self.script = "\n".join(script_chunks)

        # attempt to validate script
        self._check_script()
        logger.debug("exiting Model")


    def initialize(
        self,
        chromosome,
        length:int=1000, #length of sim in # of generations
        mut:float=1e-7, 
        recomb:float=1e-9, 
        constants:Union[None, dict]=None,
        scripts:Union[None, list]=None,
        fileout:str="shadie.trees",
        startfile:Union[None, str]=None,
        ):
        """
        Initialize a simulation. This fills the SLIM intialize() code
        block with MutationType, ElementType, Element, and other init
        type code.
        """
        logger.debug("initializing Model")
        constants = {} if constants is None else constants
        scripts = [] if scripts is None else scripts

        self.chromosome = chromosome
        self.length = length
        self.fileout = fileout
        self.startfile = startfile
        
        self.map['initialize'].append({
            'mutation_rate': mut,
            'recombination_rate': f"{recomb}, {int(chromosome.genome_size)}",
            'genome_size': chromosome.genome_size,
            'mutations': chromosome.to_slim_mutation_types(),
            'elements': chromosome.to_slim_element_types(),
            'chromosome': chromosome.to_slim_elements(),
            'constants': constants,
            'scripts': scripts,
        })

        if self.startfile:
            Model.early(
                self = self,
                time = 1,
                scripts = f"sim.readFromPopulationFile('{self.startfile}')",
                comment = "read starting populations from startfile"
                )
            Model.late(
                self = self,
                time = self.length,
                scripts = [
                "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n",
                f"sim.treeSeqOutput('{self.fileout}')"]
                )
        else:

            Model.late(
                self = self,
                time = self.length, 
                scripts = [
                "sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n",
                f"sim.treeSeqOutput('{self.fileout}')"],
                comment = "end of sim; save individuals, save .trees file",
            )

    def readfromfile(self,):
        """
        If a .trees file is provided, this will be the starting point
        of the simulation
        """


    def early(
        self, 
        time:Union[int, None], 
        scripts:Union[str, list], 
        comment:Union[str,None]=None,
        ):
        """
        Add event that happens before selection every generation.
        """
        self.map['early'].append({
            'time': time,
            'scripts': scripts,
            'comment': comment,
        })

    def repro(
        self, 
        population:Union[str, None], 
        scripts:Union[str, list], 
        comment:Union[str,None]=None,
        ):
        """
        Add event that happens before selection every generation.
        """
        self.map['reproduction'].append({
            'population': population,
            'scripts': scripts,
            'comment': comment,
        })

    def fitness(
        self, 
        mutation:Union[str, None], 
        scripts:Union[str, list], 
        idx:Union[str, None]=None,
        comment:Union[str,None]=None,
        ):
        """
        Add event that adjusts fitness values before fitness calc.
        """
        self.map['fitness'].append({
            'idx': idx,
            'mutation': mutation,
            'scripts': scripts,
            'comment': comment,
        })


    def survival(
        self, 
        population:Union[str, None], 
        scripts:Union[str, list], 
        idx:Union[str,None]=None,
        comment:Union[str,None]=None,        
        ):
        """
        Add event that adjusts fitness values before fitness calc.
        """
        self.map['survival'].append({
            'idx': idx,
            'population': population,
            'scripts': scripts,
            'comment': comment,
        })


    def late(
        self, 
        time:Union[int, None], 
        scripts:Union[str, list],
        idx:Union[str,None]=None,
        comment:Union[str,None]=None,
        ):
        """
        Add an event that happens after every generation (late) if 
        time is None, or only after a particular generation if time
        is an integer.
        """
        self.map['late'].append({
            'idx': idx,
            'time': time,
            'scripts': scripts,
            'comment': comment,
        })


    def custom(
        self, 
        scripts:str, 
        comment:Union[str,None]=None,
        ):
        """
        Add custom scripts outside without formatting by shadie. 
        Scripts must be Eidos-formatted. 
        """
        self.map['custom'].append({
            'scripts': scripts,
            'comment': comment,
        })


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
        
        #newbuild cmd 
        cmdnb = ["/Users/elissa/bin/slim", "/tmp/slim.slim"]

        # capture stdout
        proc = subprocess.Popen(cmdnb, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
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
        
        # print(chrom.data.iloc[:, :4,])
        # print(chrom.mutations)

        # init the model
        model.initialize(chromosome=chrom)
        
        model.early(
            time=1000, 
            scripts="sim.addSubpop('p1', 1000)", 
            comment="diploid sporophytes",
        )
        model.fitness(
            mutation="m4",
            scripts="return 1 + mut.selectionCoeff",
            comment="gametophytes have no dominance effects, s1",
        )

        model.custom(
            scripts="s2 fitness(m5) {\n    return 1 + mut.selectionCoeff;\n}",
            comment="gametophytes have no dominance effects",
        )


    print(model.script)
    #model.run()
