#!/usr/bin/env python

"""
A context manager for wrapping the context of a simulation setup.

Example
-------
>>> initialize() {
>>>     initializeSLiMModelType("nonWF");
>>>     defineConstant("K", 500);
>>>
>>>     initializeMutationType("m1", 0.5, "f", 0.0);
>>>     m1.convertToSubstitution = T;                      // default is T
>>>
>>>     initializeGenomicElementType("g1", m1, 1.0);
>>>     initializeGenomicElement(g1, 0, 99999);
>>>     initializeMutationRate(1e-7);
>>>     initializeRecombinationRate(1e-8);
>>> }
>>> reproduction() {
>>>     subpop.addCrossed(individual, subpop.sampleIndividuals(1));
>>> }
>>> 1 early() {
>>>     sim.addSubpop("p1", 10);
>>> }
>>> early() {
>>>     p1.fitnessScaling = K / p1.individualCount;
>>> }
>>> late() {
>>>     inds = p1.individuals;
>>>     catn(sim.generation + ": " + size(inds) + " (" + max(inds.age) + ")");
>>> }
>>> 2000 late() {
>>>     sim.outputFull(ages=T);
>>> }
"""

from typing import Union, Optional
import os
import tempfile
import subprocess
from contextlib import AbstractContextManager
from concurrent.futures import ProcessPoolExecutor
from loguru import logger
import numpy as np
from shadie.base.mutations import MutationTypeBase
from shadie.base.elements import ElementType
from shadie.reproduction.optimized.api_opt import ReproductionApi
from shadie.sims.format import (
    format_event_dicts_to_strings,
    EVENT_TO_FORMATTER,
)

# cannot do both mutationRate and nucleotidebased


class Model(AbstractContextManager):
    """The core shadie Class object to create and run SLiM scripts.

    A :class:`shadie.Model` class acts as a context manager using the
    `with` statement to wrap code to define a model, ensuring that all
    arguments are checked for compatibility and that the minimal
    requirements are present (e.g., an initialize call).

    Attributes
    ----------
    map: dict
        This dictionary stores SLiM command snippets until the model
        context closes, at which time they are converted to a string
        and stored in :attr:`.script`.
    script: str
        The final built SLiM script composing all of the arguments
        passed to the Model object. This script can be run from the
        Model object by calling :meth:`~shadie.Model.run`.
    stdout: str
        The stdout from running the SLiM script in subprocess is
        stored. If an error occurs an Exception will be raised.
    ...

    Examples
    --------
    >>> model = shadie.Model()
    >>> chrom = shadie.chromosome.default()
    >>> with shadie.Model() as model:
    >>>     model.initialize(chromosome=chrom, sim_time=10000)
    >>>     model.reproduction.bryophyte("dio", 10000, 10000)
    >>> model.run()
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
            'custom': [],
        }
        """A dict to store SLiM script snippets until context closure."""
        self.script = ""
        """The final SLiM script built from components in `.map`."""
        self.stdout = ""
        """The stdout from running `slim script` if `.run()` is called."""

        # store chromosome (mut, ele), constants, and populations
        self.sim_time: int = 0
        """The length of the simulation in sprorophyte generations."""
        self.chromosome = None
        """Chromosome object with genome structure."""
        
        self.constants = {}
        self.populations = {}

        self.reproduction = ReproductionApi(self)
        """API to access reproduction functions."""

    def __repr__(self):
        return "<shadie.Model generations={self.sim_time}>"

    def __enter__(self):
        """
        On entry the Class counters of mutation and element types
        is reset to zero. The .map dictionary is also cleared.
        """
        MutationTypeBase.idx = 0
        ElementType.idx = 0
        self.map = {i: [] for i in self.map}
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

                # check event key is valid
                if key not in EVENT_TO_FORMATTER:
                    raise ValueError(f"'{key}' not supported")

                formatter = EVENT_TO_FORMATTER[key]
                script = formatter.format(**event)
                script_chunks.append(script)

        # collapse into the final script string
        self.script = "\n".join(script_chunks)

        # attempt to validate script
        self._check_script()
        logger.debug("exiting Model")


    def initialize(
        self,
        chromosome,
        sim_time:int=1000, #length of sim in # of diploid generations
        mutation_rate:float=1e-8, #mutation rate
        recomb_rate:float=1e-9, #recombination rate
        ne:Union[None, int]=None, #option ne parameter in initialize
        file_in:Union[None, str]=None, #optional starting file
        constants:Union[None, dict]=None,
        scripts:Union[None, list]=None,
        file_out: str="shadie.trees",
        ):
        """Add an initialize() block to the SLiM code map.

        This sets the chromosome structure by initializing MutationType
        and ElementType objects, the length of the simulation in number
        of generations, sets mutation and recombination rates, and can
        add other constants or scripts, and the file paths for i/o.

        Parameters
        ----------
        chromosome: shadie.Chromosome
            A Chromosome object describes the Elements that compose
            the genome, and thus where.
        sim_time: int
            The number of *sporophyte* generations (len of simulation).
        mutation_rate: float
            The per-site per-*sporophyte*-generation mutation rate.
            This is applied at stage.
        recomb_rate: float
            The per-site per-*sporophyte*-generation recombination rate.
            This is applied in the sporophyte generation during meiosis.
        ne: int
            Optional, for WF (base) model ONLY - user can define ne here
            if they are not going to call .reproduction
        file_in: str
            Optional starting .trees file used to initialize the starting
            population
        constants: dict[str,Any]
            Custom constants defined by user
        scripts: list[str]
            Customo scripts provided by the user
        file_out: str
            Filepath to save output
        ...
        """
        logger.debug("initializing Model")
        constants = {} if constants is None else constants
        scripts = [] if scripts is None else scripts

        self.chromosome = chromosome
        self.sim_time = sim_time
        self.file_in = file_in
        self.file_out = file_out
        self.ne = ne
        self.mutation_rate = mutation_rate
        
        self.map['initialize'].append({
            'mutation_rate': mutation_rate,
            'recombination_rate': f"{recomb_rate}, {int(chromosome.genome_size)}",
            'genome_size': chromosome.genome_size,
            'mutations': chromosome.to_slim_mutation_types(),
            'elements': chromosome.to_slim_element_types(),
            'chromosome': chromosome.to_slim_elements(),
            'constants': constants,
            'scripts': scripts,
        })

        # Standard single population simulation.
    
    def readfromfile(self, tag_scripts:str):
        """
        If a .trees file is provided, this will be the starting point
        of the simulation
        """
        if self.file_in:
        # starting from another simulation starting point.
            #raise NotImplementedError("This isn't ready to use yet.")
            scripts = [f"sim.readFromPopulationFile('{self.file_in}')"]
            for i in tag_scripts:
                scripts.append(i)

            self.early(
                time=1,
                scripts=scripts,
                comment="read starting populations from file_in"
                )
            # self.late(
            #     time=1,
            #     scripts=["sim.treeSeqRememberIndividuals(sim.subpopulations.individuals)\n"],
            #     comment="save starting pop")

    def early(
        self,
        time: Union[int, None],
        scripts: Union[str, list],
        idx: Union[str,None]= None,
        comment: Union[str,None]=None,
        ):
        """Add an early() block to the SLiM code map.

        Events in `early` blocks occur before selection in every
        generation if time=None, or only in a specified generation if
        a time arg is entered.
        """
        self.map['early'].append({
            'time': time,
            'scripts': scripts,
            'idx': idx,
            'comment': comment,
        })

    def repro(
        self,
        population: Union[str, None],
        scripts: Union[str, list],
        comment: Union[str,None]=None,
        ):
        """Add a custom reproduction() block to the SLiM code map.

        Users will usually want to use the organism specific functions
        in the :meth:`Model.reproduction` API rather than write custom
        functions, since a proper reproduction cycle often involves
        entering arguments to multiple event types.
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
        time: Union[int, None],
        scripts: Union[str, list],
        idx: Union[str,None]=None,
        comment: Union[str,None]=None,
        ):
        """Add a late() block to the SLiM code map.

        Events in `late` blocks occur at the end of *every* generation
        if time=None, or only a specified generation if a time is int.
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

    def _check_repro(self):
        """
        TODO.
        """
        #check for reproduction; if does not exist, implement WF model
        if "reproduction" in self.script:
            pass
        else:
            print(self.ne)
            self.reproduction.base()
            logger.warning("You did not specify a reproduction mode "
                "so a default WF model has been used for this simulation")


    def _check_script(self):
        """
        TODO.
        """
        #check for initialization
        assert "initialize" in self.script, (
            "You must call initialize() from within Model context")

        # assert "reproduction" in self.script, (
        #     "You must specify a reproduction model to use")

        # assert "late" in self.script, (
            # "You must call late() from within Model context")

    def write(self, path: Optional[str]=None):
        """Write SLIM script to the outname filepath or stdout."""
        if path is None:
            print(self.script)
        else:
            with open(path, 'w') as out:
                out.write(self.script)

    def run(self, seed: int=None, binary: Optional[str]=None):
        """Run `slim script` using subprocess and store STDOUT to self.stdout.

        The script is written to a tmpfile and run in subprocess using
        the `slim` command in the user's PATH, or the path to a binary
        passed as an argument.

        Parameters
        ----------
        seed: Optional[int]
            A seed for the random number generator.
        binary: Optional[str]
            The full path to a slim binary.

        Examples
        --------
        >>> model = shadie.Model()
        >>> chrom = shadie.chromosome.default()
        >>> with model:
        >>>     model.initialize(chromosome=chrom, sim_time=1000)
        >>> model.run()
        """
        # get a slim binary and random number generator
        binary = binary if binary is not None else "slim"

        # find and write to a tmpdir, run it, collect result.
        with tempfile.TemporaryDirectory() as tmpdirname:

            # write slim script to tempfile
            script_path = os.path.join(tmpdirname, f"slim-{seed}.slim")
            self.write(script_path)

            # build command and run in subprocess
            cmd = [binary, '-seed', str(seed), script_path]
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # capture stdout
            out, _ = proc.communicate()

            # check for errors
            if proc.returncode:
                logger.error(out.decode())
                raise SyntaxError("SLiM3 error")

        # todo: parse stdout to store, and send warnings to logger
        self.stdout = out.decode()
        # _, warnings, stdout = out.decode().split("#")
        # self.stdout = stdout
        # logger.warning(warnings.split("//")[0].strip())
        # ...

        # TODO: idea: we could return it as a TreeSequence? or even
        # as the fully recapitated TreeSequence, though that is not
        # useful for when we want to combine two sims....


    def run_parallel(
        self,
        njobs: int=2,
        seed: Optional[int]=None,
        binary: Optional[str]=None,
        ):
        """Submit jobs to run in parallel. NOT READY YET.

        TODO: this needs to edit the trees file path to be different
        for the 2-N runs.
        """
        rng = np.random.default_rng(seed)

        rasyncs = []
        with ProcessPoolExecutor(njobs) as pool:
            for _ in range(njobs):
                args = (rng.integers(2**31), binary)
                rasync = pool.submit(self.run, *args)
                rasyncs.append(rasync)



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
