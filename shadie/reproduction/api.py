#!/usr/bin/env python

"""
API to make reproduction functions accessible from Model object.
These are the main user-facing options for implementing life 
histories into SLiM scripts using the shadie Model context.
"""

from shadie.reproduction.base import Base
from shadie.reproduction.bryobase import Bryophyte
from shadie.reproduction.angiobase import Spermatophyte
from shadie.reproduction.fernbase import Pteridophyte
from typing import Union


class ReproductionApi:
    """API for generating organism specific reproduction code blocks.

    Methods
    -------
    bryophyte
    spermatphyte
    angiosperm
    """
    def __init__(self, model: 'shadie.Model'):
        self.model = model

    def bryophyte(
        self, 
        mode:str, 
        spo_ne: int,
        gam_ne: int,
        spo_mutation_rate: Union[None, float]=None,
        gam_mutation_rate: Union[None, float]=None,
        spores_per_spo: int=100,
        gam_female_to_male_ratio: float.as_integer_ratio = (1,1),
        gam_clone_rate: float=0.0,
        gam_clone_number: int = 1,
        gam_self_rate: float=0.,
        gam_maternal_effect: float=0,
        spo_random_death_chance:float=0,
        gam_random_death_chance: float=0, 
        _file_in = None,
        _chromosome=None,
        _sim_time = None,  
        _file_out = None,          
        ):
        """Adds bryo life history to the model scripts dict.

        Generate scripts appropriate for a bryophyte (moss, liverwort,
        or hornwort) life history. This adds code to the following 
        SLiM script blocks: reproduction, early, ...

        Parameters
        ----------
        mode: str
            A life history strategy: "dio" or "mono"-icous. 
        spo_ne: int
            Sporophyte (diploid) effective population size.
        gam_ne: int
            Gametophyte (haploid) effective population size.
        spo_mutation_rate: float
            Sporophyte mutation rate; chance mutations will arise during
            the sporophyte generation
        gam_mutation_rate: float
            Gametophyte mutation rate; chance mutations will arise during
            the gametophyte generation
        spores_per_spo: int
            Number of spores generated by each spororphyte.
        gam_female_to_male_ratio: tuple
            Gametophyte female:male ratio; e.g. (1,1)
        gam_clone_rate: float
            Chance a gametophyte will clone
        gam_selfing_rate: float
            Chance a gametophyte will engage in gametophytic selfing
        gam_maternal_effect_weight: flat
            Maternal contribution to diploid offspring fitness (as
            weighted average)
        spo_random_death_chance:float
            Random chance a sporophyte will die before reproducing, 
            regardless of fitness. 
        gam_random_death_chance: float
            Random chance a gametophyte will die before reproducing,
            regardless of fitness
        file_in: str
            Provide a .trees file that will serve as the starting 
            point for the simulation
        ...
        """
        Bryophyte(
            model=self.model,
            mode=mode, spo_ne=spo_ne, gam_ne=gam_ne,
            spo_mutation_rate=spo_mutation_rate,
            gam_mutation_rate=gam_mutation_rate, 
            gam_female_to_male_ratio=gam_female_to_male_ratio,
            spores_per_spo=spores_per_spo, 
            gam_clone_rate=gam_clone_rate,
            gam_clone_number=gam_clone_number,
            gam_self_rate=gam_self_rate,
            gam_maternal_effect=gam_maternal_effect,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            _file_in=self.model.file_in, _chromosome=self.model.chromosome, 
            _sim_time = 2*self.model.sim_time, 
            _file_out = self.model.file_out,
        ).run()


    def pteridophyte(
        self, 
        mode:str, 
        spo_ne: int,
        gam_ne: int,
        spo_mutation_rate: Union[None, float]=None,
        gam_mutation_rate: Union[None, float]=None,
        spores_per_spo: int=100,
        spo_female_to_male_ratio: float.as_integer_ratio = (1,1),
        gam_female_to_male_ratio: float.as_integer_ratio = (1,1),
        spo_clone_rate: float=0.0,
        spo_clone_number: int=1,
        gam_clone_rate: float=0.0,
        gam_clone_number: int=1,
        gam_self_rate: float=0.0,
        gam_maternal_effect: float=0,
        spo_random_death_chance: float=0, 
        gam_random_death_chance: float=0,
        _file_in = None,
        _chromosome=None,
        _sim_time = None,  
        _file_out = None,       
        ):
        """
        Generate scripts appropriate for an pteridophyte (lycophytes 
        and ferns) life history.

        Parameters:
        -----------
        mode: str
            A life history strategy or "homo" or "hetero" -sporous.
        spo_ne: int
            Sporophyte (diploid) effective population size.
        gam_ne: int
            Gametophyte (haploid) effective population size.
        spo_mutation_rate: float
            Sporophyte mutation rate; chance mutations will arise during
            the sporophyte generation
        gam_mutation_rate: float
            Gametophyte mutation rate; chance mutations will arise during
            the gametophyte generation
        spores_per_spo: int
            Number of spores generated by each spororphyte.
        spo_female_to_male_ratio: tuple
            Sporophyte female:male ratio; e.g. (1,1) 
        gam_female_to_male_ratio: tuple
            Gametophyte female:male ratio; e.g. (1,1)
        spo_clone_rate: float
            Chance a sporophyte will clone
        gam_clone_rate: float
            Chance a gametophyte will clone
        gam_selfing_rate: float
            Chance a gametophyte will engage in gametophytic selfing
        gam_maternal_effect_weight: flat
            Maternal contribution to diploid offspring fitness (as
            weighted average)
        spo_random_death_chance:float
            Random chance a sporophyte will die before reproducing, 
            regardless of fitness. 
        gam_random_death_chance: float
            Random chance a gametophyte will die before reproducing,
            regardless of fitness
        file_in: str
            Provide a .trees file that will serve as the starting 
            point for the simulation
        ...
        """
        Pteridophyte(
            model=self.model, mode=mode, spo_ne=spo_ne, gam_ne=gam_ne,
            spo_mutation_rate = spo_mutation_rate,
            gam_mutation_rate = gam_mutation_rate, 
            spo_female_to_male_ratio = spo_female_to_male_ratio,
            gam_female_to_male_ratio = gam_female_to_male_ratio,
            spores_per_spo=spores_per_spo,
            spo_clone_rate=spo_clone_rate, 
            spo_clone_number = spo_clone_number,
            gam_clone_rate=gam_clone_rate,
            gam_clone_number = gam_clone_number,
            gam_self_rate=gam_self_rate,
            gam_maternal_effect=gam_maternal_effect,
            spo_random_death_chance=spo_random_death_chance, 
            gam_random_death_chance=gam_random_death_chance,
           _file_in=self.model.file_in,
            _chromosome=self.model.chromosome,
            _sim_time = 2*self.model.sim_time, 
            _file_out = self.model.file_out,
        ).run()


    def spermatophyte(
        self, 
        mode:str, 
        spo_ne: int,
        gam_ne: int,
        spo_mutation_rate: Union[None, float]=None,
        gam_mutation_rate: Union[None, float]=None,
        spo_female_to_male_ratio: float.as_integer_ratio = (1,1),
        spo_clone_rate: float=0.0,
        spo_clone_number:  int =1,
        ovule_count: int=30,
        fertilization_rate: float=0.7,
        pollen_count: int=100,
        pollen_comp: str="F",
        pollen_per_stigma: int=5,
        spo_random_death_chance: float=0,  
        gam_random_death_chance: float=0, 
        _file_in = None,
        _chromosome=None,
        _sim_time = None,  
        _file_out = None, 
        ):
        """
        Generate scripts appropriate for an angiosperm (flowering plant)
        life history. Appropriate for gymnosperms as well.
        This adds code to the following 
        SLiM script blocks: reproduction, early, ...

        Parameters:
        -----------
        mode: str
            A life history strategy or "dio" or "mono" -ecious.
        ...
        """
        Spermatophyte(
            model=self.model, mode=mode,spo_ne=spo_ne, gam_ne=gam_ne,
            spo_mutation_rate = spo_mutation_rate,
            gam_mutation_rate = gam_mutation_rate, 
            spo_female_to_male_ratio = spo_female_to_male_ratio,
            spo_clone_rate=spo_clone_rate, 
            spo_clone_number = spo_clone_number, ovule_count=ovule_count,
            fertilization_rate=fertilization_rate, 
            pollen_count=pollen_count, pollen_comp=pollen_comp,
            pollen_per_stigma=pollen_per_stigma,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            _file_in=self.model.file_in,
            _chromosome=self.model.chromosome, 
            _sim_time = 2*self.model.sim_time, 
            _file_out = self.model.file_out,
        ).run()

    def base(
        self,  
        ne = None,
        sexes = False,    
        ):
        """
        Generate scripts appropriate for basic SLiM nonWF model, set up
        as a WF model

        Parameters:
        -----------
        sexes: bool
            default = False; individuals are hemraphroditic. If True, 
            individuals will be male and female
        ...
        """
        Base(
            model=self.model, ne = ne, sexes = sexes, 
            _chromosome=self.model.chromosome, 
            _sim_time = self.model.sim_time, 
            _file_out = self.model.file_out,
        ).run()

