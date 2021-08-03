#!/usr/bin/env python

"""
API to make reproduction functions accessible from Model object.
These are the main user-facing options for implementing life 
histories into SLiM scripts using the shadie Model context.
"""

from shadie.reproduction.base import Bryophyte
from shadie.reproduction.angiobase import Angiosperm
from shadie.reproduction.fernbase import Pteridophyte



class ReproductionApi:
    """
    Reproduction API for accessing functions to generate organism
    specific reproduction code blocks.

    - model.reproduction.bryophyte()
    - model.reproduction.spermatphyte()
    - model.reproduction....()

    """
    def __init__(self, model: 'shadie.Model'):
        self.model = model


    def bryophyte(
        self, 
        chromosome,
        mode:str, 
        diploid_ne: int,
        haploid_ne: int,
        female_to_male_ratio: float=0.5,
        spores_per_sporophyte: int=100,
        clone_rate: float=0.0,
        selfing_rate: float=0.,
        maternal_effect_weight: float=0,
        random_death_chance: float=0,           
        ):
        """
        Generate scripts appropriate for a bryophyte (moss, liverwort,
        or hornwort) life history. This adds code to the following 
        SLiM script blocks: reproduction, early, ...

        Parameters:
        -----------
        mode: str
            A life history strategy or "dio" or "mono" -icous.
        ...
        """
        Bryophyte(
            model=self.model, chromosome=chromosome, mode=mode,
            diploid_ne=diploid_ne, haploid_ne=haploid_ne,
            female_to_male_ratio=female_to_male_ratio,
            spores_per_sporophyte=spores_per_sporophyte,
            clone_rate=clone_rate,
            selfing_rate=selfing_rate,
            maternal_effect_weight=maternal_effect_weight,
            random_death_chance=random_death_chance,
        ).run()


    def pteridophyte(
        self, 
        chromosome,
        mode:str, 
        diploid_ne: int,
        haploid_ne: int,
        spores_per_sporophyte: int=100,
        female_to_male_ratio: float=0.5,
        gam_female_to_male_ratio: float=0.5,
        clone_rate: float=0.0,
        gam_clone_rate: float=0.0,
        selfing_rate: float=0.0,
        maternal_effect_weight: float=0,
        random_death_chance: float=0,       
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
        Pteridophyte(
            model=self.model, chromosome=chromosome, mode=mode,
            diploid_ne=diploid_ne, haploid_ne=haploid_ne,
            female_to_male_ratio=female_to_male_ratio,
            gam_female_to_male_ratio = gam_female_to_male_ratio,
            spores_per_sporophyte=spores_per_sporophyte,
            clone_rate=clone_rate, gam_clone_rate=gam_clone_rate,
            selfing_rate=selfing_rate,
            maternal_effect_weight=maternal_effect_weight,
            random_death_chance=random_death_chance,
        ).run()

    def spermatophyte(self, ):
        """
        TODO:
        """

    def angiosperm(
        self, 
        chromosome,
        mode:str, 
        diploid_ne: int,
        haploid_ne: int,
        female_to_male_ratio: float=0.5,
        clone_rate: float=0.0,
        ovule_count: int=30,
        fertilization_rate: float=0.7,
        pollen_count: int=100,
        pollen_comp: str="F",
        pollen_per_stigma: int=5,
        random_death_chance: float=0,       
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
        Angiosperm(
            model=self.model, chromosome=chromosome, mode=mode,
            diploid_ne=diploid_ne, haploid_ne=haploid_ne,
            female_to_male_ratio=female_to_male_ratio,
            clone_rate=clone_rate, ovule_count=ovule_count,
            fertilization_rate=fertilization_rate, 
            pollen_count=pollen_count, pollen_comp=pollen_comp,
            pollen_per_stigma=pollen_per_stigma,
            random_death_chance=random_death_chance,
        ).run()

