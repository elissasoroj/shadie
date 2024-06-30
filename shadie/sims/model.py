#!/usr/bin/env python

"""A context manager for wrapping the context of a simulation setup.

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

from typing import Union, Optional, List
import os
import tempfile
import subprocess
from copy import deepcopy
from contextlib import AbstractContextManager
from concurrent.futures import ProcessPoolExecutor
from loguru import logger
import numpy as np
from shadie.base.mutations import MutationType
from shadie.base.elements import ElementType
from shadie.reproduction.api import ReproductionApi
from shadie.sims.format import (
    format_event_dicts_to_strings,
    EVENT_TO_FORMATTER,
)

# register logger to this module only
logger = logger.bind(name="shadie")


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
            'first':[],
            'reproduction': [],
            'early': [],
            'late': [],
            'muteffect': [],
            'fitness': [],
            'survival': [],
            'custom': [],
        }
        """: A dict to store SLiM script snippets until context closure."""
        self.script = ""
        """: The final SLiM script built from components in `.map`."""
        self.stdout = ""
        """: The stdout from running `slim script` if `.run()` is called."""
        self.sim_time: int = 0
        """: The length of the simulation in sporophyte generations."""
        self.chromosome = None
        """: Chromosome object with genome structure."""
        self.metadata: dict = {}
        """: Dictionary storing simulation metadata"""

        # hidden attributes set by .initialize()
        self.metadata = {
            'file_in': None,
            'file_out': None,
            'mutation_rate': None,
            'recomb_rate': None,
        }

        self.reproduction = ReproductionApi(self)
        """API to access reproduction functions."""

    def __repr__(self):
        return f"<shadie.Model generations={self.sim_time}>"

    def __enter__(self):
        """
        On entry the Class counters of mutation and element types
        is reset to zero. The .map dictionary is also cleared.
        """
        MutationType.idx = 0
        ElementType.idx = 0
        self.map = {i: [] for i in self.map}
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """
        Fills the script with context-defined content in a set
        order and runs checks on the script.
        """
        sorted_keys = [
            'initialize', 'shadie', 'first', 'reproduction', 'early',
            'custom', 'muteffect', 'survival', 'fitness', 'late',
        ]

        # copy map and split timed events to a new key list
        mapped = self.map.copy()
        mapped['shadie'] = []
        mapped_keys = list(mapped.keys())
        for key in mapped_keys:
            if key == "custom":
                for i in range(0, len(mapped[key])):
                    script = mapped[key][i]
                    if "DEFINITIONS" in script['comment']:
                        mapped['shadie'].append(mapped[key].pop(i))

        #     for item in mapped[key]:
        #         print(item)
        #         if 'time' in item:
        #             if item['time']==1:
        #                 mapped['1'].append(mapped[key].pop('time' == 1))
        #             else:
        #                 mapped['timed'].append(mapped[key].pop(0))

        # logger.info(f"{test}")
        # visit events by ordered key type
        script_chunks = []
        for key in sorted_keys:
            # sort events within key type
            if key == "shadie":
                events = mapped[key]
            elif key == "first":
                events = sorted(mapped[key], key=lambda x: str(x['time']))
            elif key == "early":
                events = sorted(mapped[key], key=lambda x: str(x['time']))
            elif key == "late":
                events = sorted(mapped[key], key=lambda x: str(x['time']))
            elif key == "muteffect":
                events = sorted(mapped[key], key=lambda x: -1 * float('inf') if x['idx'] is None else x['idx'])
            elif key == "fitness":
                events = sorted(mapped[key], key=lambda x: -1 * float('inf') if x['idx'] is None else x['idx'])
            elif key == "survival":
                events = sorted(mapped[key], key=lambda x: -1 * float('inf') if x['idx'] is None else x['idx'])
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
        sim_time: int = 1000,
        mutation_rate: float = 1e-8,
        recomb_rate: float = 1e-9,
        constants: Union[None, dict] = None,
        simglobals: Union[None, dict] = None,
        scripts: Union[None, list] = None,
        file_in: Union[None, str] = None,
        file_out: str = "shadie.trees",
        skip_neutral_mutations: bool = False,
        _simplification_interval: Union[str, int] = "NULL",
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
        constants: dict[str,Any]
            Custom constants defined by user
        simglobals: dict[str,Any]
            Custom globals defined by user
        scripts: list[str]
            Customo scripts provided by the user
        file_in: str
            Optional starting .trees file used to initialize the
            starting population
        file_out: str
            Filepath to save output
        skip_neutral_mutations: bool
            If True then mutations are not added to neutral genomic
            regions. This should be used if you plan to add coalescent
            recapitated ancestry and mutations. Default=False.
        """
        logger.debug("initializing Model")
        constants = {} if constants is None else constants
        simglobals = {} if simglobals is None else simglobals
        scripts = [] if scripts is None else scripts

        # store a copy of the chromosome and set to keep or exclude neutral.
        self.chromosome = deepcopy(chromosome)
        self.chromosome._skip_neutral_mutations = skip_neutral_mutations
        self.sim_time = sim_time
        self.metadata.update({
            'file_in': file_in,
            'file_out': file_out,
            'mutation_rate': mutation_rate,
            'recomb_rate': recomb_rate,
        })

        self.map['initialize'].append({
            'simplification_interval': _simplification_interval,
            'mutation_rate': mutation_rate,
            'recombination_rate': f"{recomb_rate}, {int(self.chromosome.genome_size)}",
            'genome_size': self.chromosome.genome_size,
            'mutations': self.chromosome.to_slim_mutation_types(),
            'mutation_names': self.chromosome.mutation_list(),
            'elements': self.chromosome.to_slim_element_types(),
            'chromosome': self.chromosome.to_slim_elements(),
            'constants': constants,
            'simglobals': simglobals,
            'scripts': scripts,
        })

    def _read_from_file(self, tag_scripts: List[str]):
        """Set an existing .trees file as starting state of simulation.

        If the trees file is not a shadie trees file (e.g., with
        subpops defined as p0 and p1) this will cause problems.
        """
        scripts = [f"sim.readFromPopulationFile('{self.metadata['file_in']}')"]
        scripts.extend(tag_scripts)
        self.early(
            time=1,
            scripts=scripts,
            comment="read starting populations from file_in"
        )

    def first(
        self,
        time: Union[int, None],
        scripts: Union[str, list],
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
    ):
        """Add an first() block to the SLiM code map.

        Events in `first` blocks occur before reproduction in every
        generation if time=None, or only in a specified generation if
        a time arg is entered.
        """
        logger.debug(f"define first() @time={time} for idx={idx}")
        self.map['first'].append({
            'time': time,
            'scripts': scripts,
            'idx': idx,
            'comment': comment,
        })

    def early(
        self,
        time: Union[int, None],
        scripts: Union[str, list],
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
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
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
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
            'idx': idx,
            'comment': comment,
        })

    def muteffect(
        self,
        mutation: Union[str, None],
        scripts: Union[str, list],
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
    ):
        """Add event that adjusts fitness values before fitness calc.

        Parameters
        ----------
        mutation: str or None
            The name of a MutationType object (e.g., MutationType.name)
        scripts: str or list
            One or more SLiM scripts to execute.
        idx: str or None
            An index...
        comment: str or None
            An optional comment to embed in SLiM code.

        Example
        -------
        >>> model.muteffect(
        >>>     idx=None, mutation=mut.name, scripts=SCRIPT,
        >>>     comment="mutation only affects haploid",
        >>> )
        """
        logger.debug(f"define mutationEffect() for mutationType={mutation}")
        self.map['muteffect'].append({
            'idx': idx,
            'mutation': mutation,
            'scripts': scripts,
            'comment': comment,
        })

    def fitness(
        self,
        population: str,  # subpop or ind
        scripts: Union[str, list],
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
    ):
        """Add event that adjusts fitness values before fitness calc.
        """
        self.map['fitness'].append({
            'idx': idx,
            'population': population,
            'scripts': scripts,
            'comment': comment,
        })

    def survival(
        self,
        population: Union[str, None],
        scripts: Union[str, list],
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
    ):
        """Add event that adjusts fitness values before fitness calc.
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
        idx: Union[str, None] = None,
        comment: Union[str, None] = None,
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
        scripts: str,
        comment: Union[str, None] = None,
    ):
        """Add custom scripts outside without formatting by shadie.
        Scripts must be Eidos-formatted.
        """
        self.map['custom'].append({
            'scripts': scripts,
            'comment': comment,
        })

    def _check_script(self):
        """Checks that the model contains the minimal requirements.
        """
        assert "initialize" in self.script, (
            "You must call initialize() from within Model context.")
        assert "reproduction" in self.script, (
            "You must call reproduction() from within Model context "
            "to implement either an organism specific reproduction "
            "or a standard wright_fisher model.")

    def write(self, path: Optional[str] = None) -> None:
        """Write SLIM script to the outname filepath or stdout."""
        if path is None:
            print(self.script)
        else:
            with open(path, 'w') as out:
                out.write(self.script)

    def run(self, seed: int = None, binary: Optional[str] = None):
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
            with subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            ) as proc:

                # capture stdout
                out, _ = proc.communicate()

                # check for errors
                if proc.returncode:
                    logger.error(out.decode())
                    self.write("/tmp/slim.slim")
                    raise SyntaxError("SLiM error, see script at /tmp/slim.slim")

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
        njobs: int = 2,
        seed: Optional[int] = None,
        binary: Optional[str] = None,
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
    shadie.set_log_level = "DEBUG"

    with shadie.Model() as model:

        # define mutation types
        m0 = shadie.mtype(0.5, 'n', [2.0, 1.0])
        m1 = shadie.mtype(0.5, 'g', [3.0, 1.0], affects_diploid=False)

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
        print(chrom.mutations)

        # init the model
        model.initialize(chromosome=chrom, file_out="/tmp/shadie.trees")
        model.reproduction.bryophyte_monoicous(
            spo_pop_size=1000,
            gam_pop_size=2000,
        )

        model.early(
            time=1000,
            scripts="sim.addSubpop('p1', 1000)",
            comment="diploid sporophytes",
        )

        model.first(
            time=100,
            scripts="sim.addSubpop('p0', 100)",
            comment="add haploid gametos",
        )

        model.fitness(
            population='p1',
            scripts="return 1 + mut.selectionCoeff",
            comment="gametophytes have no dominance effects",
            idx=1
        )

        # model.custom(
        #     scripts="s2 fitness(m5) {\n    return 1 + mut.selectionCoeff;\n}",
        #     comment="gametophytes have no dominance effects",
        # )

    # print(model.script)
    model.run()
