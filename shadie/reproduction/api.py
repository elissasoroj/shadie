#!/usr/bin/env python

"""
API to make reproduction functions accessible from Model object
"""

from shadie.reproduction.base import Bryophyte



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
        clone_rate: float=1.0,
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
        """
        Bryophyte(
            model=self.model, chromosome = chromosome, mode=mode,
            diploid_ne=diploid_ne, haploid_ne=haploid_ne,
            female_to_male_ratio=female_to_male_ratio,
            spores_per_sporophyte=spores_per_sporophyte,
            clone_rate=clone_rate,
            selfing_rate=selfing_rate,
            maternal_effect_weight=maternal_effect_weight,
            random_death_chance=random_death_chance,
        ).run()
