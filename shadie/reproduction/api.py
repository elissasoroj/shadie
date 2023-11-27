#!/usr/bin/env python

"""
API to make reproduction functions accessible from Model object.
These are the main user-facing options for implementing life
histories into SLiM scripts using the shadie Model context.
"""

from typing import Union, Tuple, Optional

from shadie.reproduction.base import WrightFisher
from shadie.reproduction.wfspecial import ClonalHaploidWF, AltGenWF, HaploidWF
from shadie.reproduction.bryobase import BryophyteMonoicous, BryophyteDioicous
from shadie.reproduction.angiobase import AngiospermMonoecious, AngiospermDioecious
from shadie.reproduction.fernbase import PteridophyteHomosporous, PteridophyteHeterosporous


class ReproductionApi:
    """API for generating organism specific reproduction code blocks.

    This is called from an initialized Model object to access functions
    as Model.reproduction.bryophyte(...), for example.
    """
    def __init__(self, model: 'shadie.Model'):
        self.model = model

    def bryophyte_monoicous(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        spo_mutation_rate: Optional[float]=None,
        gam_mutation_rate: Optional[float]=None,
        gam_clone_rate: float=0.8,
        gam_clones_per: int=100,
        #gam_self_rate: float=0.2,
        gam_self_rate_per_egg: float=0.2,
        spo_self_rate_per_egg: float=0.0,
        #spo_self_rate: float=0.0,
        gam_maternal_effect: float=0.0,
        spo_random_death_chance: float=0.0,
        gam_random_death_chance: float=0.0,
        spo_spores_per: int=500,
        gam_archegonia_per: int=5,
        gam_k: int=None,
        ):
        """Adds Monoicous Bryophyte life history to model.

        Generate scripts appropriate for a bryophyte (moss, liverwort,
        or hornwort) life history. This adds code to the following
        SLiM script blocks: reproduction, early, ...

        Parameters
        ----------
        spo_pop_size: int
            Sporophyte (diploid) population carrying capacity.
        gam_pop_size: int
            Gametophyte (haploid) population carrying capacity.
        spo_mutation_rate: float
            Sporophyte mutation rate; mutations per site per generation
            applied each sporophyte generation. If rates are set for 
            spo and gam separately here then each is set to 1/2 the 
            initialize mutation rate. 
        gam_mutation_rate: float
            Gametophyte mutation rate; mutations per site per generation
            applied each gametophyte generation.
        spo_spores_per: int
            Number of spores generated by each sporophyte.
        gam_clone_rate: float
            Chance a gametophyte will clone
        gam_selfing_rate: float
            Chance a gametophyte will engage in gametophytic selfing
        gam_maternal_effect_weight: float
            Maternal contribution to diploid offspring fitness (as
            weighted average)
        spo_random_death_chance:float
            Random chance a sporophyte will die before reproducing,
            regardless of fitness.
        gam_random_death_chance: float
            Random chance a gametophyte will die before reproducing,
            regardless of fitness
        ...
        """
        BryophyteMonoicous(
            model=self.model,
            spo_pop_size=spo_pop_size,
            gam_pop_size=gam_pop_size,
            spo_mutation_rate=spo_mutation_rate,
            gam_mutation_rate=gam_mutation_rate,
            gam_clone_rate=gam_clone_rate,
            gam_clones_per=gam_clones_per,
            #gam_self_rate=gam_self_rate,
            gam_self_rate_per_egg=gam_self_rate_per_egg,
            spo_self_rate_per_egg=spo_self_rate_per_egg,
            #spo_self_rate=spo_self_rate,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            gam_maternal_effect=gam_maternal_effect,
            spo_spores_per=spo_spores_per,
            gam_archegonia_per=gam_archegonia_per,
            gam_k=gam_k,
        ).run()

    def bryophyte_dioicous(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        spo_mutation_rate: Optional[float]=None,
        gam_mutation_rate: Optional[float]=None,
        gam_female_to_male_ratio: Tuple[float,float]=(2, 1),
        gam_clone_rate: float=0.8,
        gam_clones_per: int=10,
        #gam_self_rate:float=0.0,
        #gam_self_rate_per_egg:float=0.0,
        spo_self_rate_per_egg: float=0.1,
        #spo_self_rate: float=0.1,
        spo_random_death_chance: float=0.08,
        gam_random_death_chance: float=0.08,
        gam_maternal_effect: float=0.5,
        spo_spores_per: int=10,
        gam_archegonia_per = 10,
        gam_k=None,
        ):
        """Adds Dioicous Bryophyte life history to model.

        Generate scripts appropriate for a bryophyte (moss, liverwort,
        or hornwort) life history. This adds code to the following
        SLiM script blocks: reproduction, early, ...

        Parameters
        ----------
        spo_pop_size: int
            Sporophyte (diploid) population carrying capacity.
        gam_pop_size: int
            Gametophyte (haploid) population carrying capacity.
        spo_mutation_rate: float
            Sporophyte mutation rate; mutations per site per generation
            applied each sporophyte generation. If rates are set for 
            spo and gam separately here then each is set to 1/2 the 
            initialize mutation rate. 
        gam_mutation_rate: float
            Gametophyte mutation rate; mutations per site per generation
            applied each gametophyte generation.
        gam_female_to_male_ratio: tuple
            Gametophyte female:male ratio; e.g. (1,1)            
        spo_spores_per: int
            Number of spores generated by each sporophyte.
        spo_megaspores_per: int
            Number of megaspores generated by each sporophyte.
        gam_sporophytes_per: int
            Number of sporophytes generated per gametophyte.
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
        ...
        """
        BryophyteDioicous(
            model=self.model,
            spo_pop_size=spo_pop_size,
            gam_pop_size=gam_pop_size,
            spo_mutation_rate=spo_mutation_rate,
            gam_mutation_rate=gam_mutation_rate,
            gam_female_to_male_ratio=gam_female_to_male_ratio,
            gam_clone_rate=gam_clone_rate,
            gam_clones_per=gam_clones_per,
            #gam_self_rate=gam_self_rate,
            #gam_self_rate_per_egg=gam_self_rate_per_egg,
            spo_self_rate_per_egg=spo_self_rate_per_egg,
            #spo_self_rate=spo_self_rate,
            gam_maternal_effect=gam_maternal_effect,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            spo_spores_per=spo_spores_per,
            gam_archegonia_per= gam_archegonia_per,
            gam_k=gam_k,
        ).run()        

    def pteridophyte_homosporous(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        spo_mutation_rate: Union[None, float]=None,
        gam_mutation_rate: Union[None, float]=None,
        gam_female_to_male_ratio: Tuple[float,float]=(2, 1),
        spo_clone_rate: float=0.0,
        spo_clones_per: int=1,
        gam_clone_rate: float=0.0,
        gam_clones_per: int=1,
        gam_self_rate_per_egg: float=0.1,
        spo_self_rate_per_egg: float=0.1,
        spo_self_rate: float=0.0,
        gam_self_rate: float=0.0,
        spo_random_death_chance: float=0,
        gam_random_death_chance: float=0,
        spo_maternal_effect: float=0,
        gam_maternal_effect: float=0,
        spo_spores_per: int=100,
        gam_archegonia_per = 10,
        gam_k: int =None,
        ):
        """
        Generate scripts appropriate for an pteridophyte (lycophytes
        and ferns) life history.

        Parameters:
        -----------
        mode: str
            A life history strategy or "homo" or "hetero" -sporous.
        spo_pop_size: int
            Sporophyte (diploid) effective population size.
        gam_pop_size: int
            Gametophyte (haploid) effective population size.
        spo_mutation_rate: float
            Sporophyte mutation rate; chance mutations will arise during
            the sporophyte generation
        gam_mutation_rate: float
            Gametophyte mutation rate; chance mutations will arise during
            the gametophyte generation
        spo_spores_per: int
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
        PteridophyteHomosporous(
            model=self.model, 
            spo_pop_size=spo_pop_size, 
            gam_pop_size=gam_pop_size,
            gam_female_to_male_ratio=gam_female_to_male_ratio,
            spo_mutation_rate=spo_mutation_rate,
            gam_mutation_rate=gam_mutation_rate,
            spo_clone_rate=spo_clone_rate,
            spo_clones_per = spo_clones_per,
            gam_clone_rate=gam_clone_rate,
            gam_clones_per=gam_clones_per,
            gam_self_rate_per_egg=gam_self_rate_per_egg,
            spo_self_rate_per_egg=spo_self_rate_per_egg,
            spo_self_rate=spo_self_rate,
            gam_self_rate=gam_self_rate,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            gam_maternal_effect=gam_maternal_effect,
            spo_maternal_effect=spo_maternal_effect,
            spo_spores_per=spo_spores_per,
            gam_archegonia_per=gam_archegonia_per,
            gam_k=gam_k,
        ).run()

    def pteridophyte_heterosporous(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        gam_female_to_male_ratio: Tuple[float,float]=(2, 1),
        spo_mutation_rate: Union[None, float]=None,
        gam_mutation_rate: Union[None, float]=None,
        spo_clone_rate: float=0.0,
        spo_clones_per: int=1,
        spo_self_rate_per_egg: float=0.1,
        spo_self_rate: float=0.0,
        spo_random_death_chance: float=0,
        gam_random_death_chance: float=0,
        spo_maternal_effect: float=0,
        gam_archegonia_per = 10,
        spo_spores_per: int=100,
        # rs_megasporangia_per: int=1,
        # rs_microsporangia_per: int=1,
        # megasporangia_megaspores_per: int=1,
        # microsporangia_microspores_per: int=100,
        gam_k: int =None,
        ):
        """
        Generate scripts appropriate for an pteridophyte (lycophytes
        and ferns) life history.

        Parameters:
        -----------
        mode: str
            A life history strategy or "homo" or "hetero" -sporous.
        spo_pop_size: int
            Sporophyte (diploid) effective population size.
        gam_pop_size: int
            Gametophyte (haploid) effective population size.
        spo_mutation_rate: float
            Sporophyte mutation rate; chance mutations will arise during
            the sporophyte generation
        gam_mutation_rate: float
            Gametophyte mutation rate; chance mutations will arise during
            the gametophyte generation
        spo_spores_per: int
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
        PteridophyteHeterosporous(
            model=self.model, 
            spo_pop_size=spo_pop_size, 
            gam_pop_size=gam_pop_size,
            gam_female_to_male_ratio=gam_female_to_male_ratio,
            spo_mutation_rate = spo_mutation_rate,
            gam_mutation_rate = gam_mutation_rate,
            spo_clone_rate=spo_clone_rate,
            spo_clones_per = spo_clones_per,
            spo_self_rate_per_egg=spo_self_rate_per_egg,
            spo_self_rate =  spo_self_rate,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            spo_maternal_effect=spo_maternal_effect,
            gam_archegonia_per=gam_archegonia_per,
            spo_spores_per=spo_spores_per,
            # rs_megasporangia_per=rs_megasporangia_per,
            # rs_microsporangia_per=rs_microsporangia_per,
            # megasporangia_megaspores_per=megasporangia_megaspores_per,
            # microsporangia_microspores_per=microsporangia_microspores_per,
            gam_k=gam_k,
        ).run()


    def angiosperm_monoecious(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        spo_mutation_rate: Optional[float]=None,
        gam_mutation_rate: Optional[float]=None,
        spo_clone_rate: float=0.0,
        spo_clones_per: int=3,
        spo_self_rate_per_egg: float=0.1,
        spo_self_rate: float=0.0,
        spo_random_death_chance: float=0.0,
        gam_random_death_chance: float=0.0,
        spo_maternal_effect: float=0.0,
        spo_pollen_per: int=100,
        spo_archegonia_per: int=30,
        fertilization_rate: float=0.7,
        pollen_success_rate: float=1.0,
        pollen_competition: str="F",
        stigma_pollen_per: int=5,
        ):
        AngiospermMonoecious(
            model=self.model, 
            spo_pop_size=spo_pop_size,
            gam_pop_size=gam_pop_size,
            spo_mutation_rate=spo_mutation_rate,
            gam_mutation_rate=gam_mutation_rate,
            spo_clone_rate=spo_clone_rate,
            spo_clones_per=spo_clones_per,
            spo_self_rate_per_egg=spo_self_rate_per_egg,
            spo_self_rate=spo_self_rate,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            spo_maternal_effect=spo_maternal_effect,
            spo_pollen_per=spo_pollen_per,
            spo_archegonia_per=spo_archegonia_per,
            fertilization_rate=fertilization_rate,
            pollen_success_rate=pollen_success_rate,
            pollen_competition=pollen_competition,
            stigma_pollen_per=stigma_pollen_per,
        ).run()
    

    def angiosperm_dioecious(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        spo_female_to_male_ratio: Tuple[float,float],
        spo_mutation_rate: Optional[float]=None,
        gam_mutation_rate: Optional[float]=0.0,
        spo_clone_rate: float=0.0,
        spo_clones_per: int=3,
        spo_random_death_chance: float=0.0,
        gam_random_death_chance: float=0.0,
        spo_maternal_effect: float=0.0,
        spo_pollen_per: int=100,
        spo_archegonia_per: int=30,
        fertilization_rate: float=0.7,
        pollen_success_rate: float=1.0,
        pollen_competition: str="F",
        stigma_pollen_per: int=5,
        ):

        AngiospermDioecious(
            model=self.model, 
            spo_pop_size=spo_pop_size,
            gam_pop_size=gam_pop_size,
            spo_female_to_male_ratio=spo_female_to_male_ratio,
            spo_mutation_rate=spo_mutation_rate,
            gam_mutation_rate=gam_mutation_rate,
            spo_clone_rate=spo_clone_rate,
            spo_clones_per=spo_clones_per,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            spo_maternal_effect=spo_maternal_effect,
            spo_pollen_per=spo_pollen_per,
            spo_archegonia_per=spo_archegonia_per,
            fertilization_rate=fertilization_rate,
            pollen_success_rate=pollen_success_rate,
            pollen_competition=pollen_competition,
            stigma_pollen_per=stigma_pollen_per,
        ).run()
    # def gymnosperm_monosporous(...)
    # def gymnosperm_heterosporous(...)    
    
    def spermatophyte(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        spo_mutation_rate: Optional[float]=None,
        gam_mutation_rate: Optional[float]=0.0,
        spo_female_to_male_ratio: float.as_integer_ratio = (1,1),
        spo_clone_rate: float=0.0,
        spo_clones_per:  int =1,
        spo_maternal_effect: float =0.0,
        spo_random_death_chance: float=0,
        gam_random_death_chance: float=0,
        spo_pollen_per: int=100,
        spo_archegonia_per: int=30,
        fertilization_rate: float=0.7,
        pollen_success_rate: float=1.0,
        pollen_competition: str="F",
        stigma_pollen_per: int=5,
        _chromosome=None,
        _sim_time = None,
        _file_in = None,
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
        spo_pop_size: int
            Sporophyte (diploid) effective population size.
        gam_pop_size: int
            Gametophyte (haploid) effective population size.
        spo_mutation_rate: float
            Sporophyte mutation rate; chance mutations will arise during
            the sporophyte generation
        gam_mutation_rate: float
            Gametophyte mutation rate; chance mutations will arise during
            the gametophyte generation. Default = 0
        spo_female_to_male_ratio: tuple
            Sporophyte female:male ratio; e.g. (1,1)
        gam_female_to_male_ratio: tuple
            Gametophyte female:male ratio; e.g. (1,1)
        spo_clone_rate: float
            Chance a sporophyte will clone
        spo_clones_per: int
            Number of clones produced by each clonal sporophyte
        ovule_count: int
            Number of ovules per sporophyte (remmeber to multiply number
            of flowers on each individual by number of ovules)
        ovule_fertilization_rate: float
            chance an ovule will be viable and set seed
        pollen_succcess_rate: float
            chance a give pollen will succcessfully fertilize an ovule
        pollen_count: int
            number of pollen produced by each sporophyte
        pollen_comp: str="F" or "T"
            turn pollen competition on or off
        pollen_per_ovule: int
            number of pollen that will compete to fertilize a single
            ovule *TODO: update code so that they compete for ALL the
            ovules in a given flower
        spo_maternal_effect_weight: float
            Maternal contribution to haploid offspring fitness (as
            weighted average)
        spo_random_death_chance:float
            Random chance a sporophyte will die before reproducing,
            regardless of fitness.
        gam_random_death_chance: float
            Random chance a gametophyte will die before reproducing,
            regardless of fitness
        ...
        """
        Spermatophyte(
            model=self.model,spo_pop_size=spo_pop_size, gam_pop_size=gam_pop_size,
            spo_mutation_rate = spo_mutation_rate,
            gam_mutation_rate = gam_mutation_rate,
            spo_female_to_male_ratio = spo_female_to_male_ratio,
            spo_clone_rate=spo_clone_rate,
            spo_clones_per = spo_clones_per, ovule_count=ovule_count,
            ovule_fertilization_rate=ovule_fertilization_rate,
            pollen_success_rate=pollen_success_rate,
            pollen_count=pollen_count, pollen_comp=pollen_comp,
            pollen_per_ovule=pollen_per_ovule,
            spo_random_death_chance=spo_random_death_chance,
            gam_random_death_chance=gam_random_death_chance,
            _chromosome=self.model.chromosome,
            _sim_time = 2*self.model.sim_time,
            _file_in=self.model.file_in,
            _file_out = self.model.file_out,
        ).run()


    def wright_fisher(
        self,
        pop_size: int,
        sexes: bool=False,  # mono/hetero terms is more consistent with others...
        ):
        """
        Generate scripts appropriate for basic WF model, set up as a 
        SLiM nonWF model

        Parameters:
        -----------
        pop_size: 
            Size of the population in number of diploids.
        sexes: bool
            default = False; individuals are hemraphroditic. If True,
            individuals will be male and female
        """
        WrightFisher(
            model=self.model,
            pop_size=pop_size,
            sexes=sexes,
        ).run()

    def wright_fisher_haploid_sexual(
        self,
        pop_size: int,
        sexes: bool=False,  # mono/hetero terms is more consistent with others...
        ):
        """
        Generate scripts appropriate for basic SLiM nonWF model, set up
        as a WF model.

        Parameters:
        -----------
        pop_size: 
            Size of the population in number of haploids.
        sexes: bool
            default = False; individuals are hemraphroditic. If True,
            individuals will be male and female
        """
        HaploidWF(
            model=self.model,
            pop_size=pop_size,
            sexes=sexes,
        ).run()

    def wright_fisher_haploid_clonal(
        self,
        pop_size: int,
        sexes: bool=False,  # mono/hetero terms is more consistent with others...
        ):
        """
        Generate scripts appropriate for basic SLiM nonWF model, set up
        as a WF model.

        Parameters:
        -----------
        pop_size: 
            Size of the population in number of haploids.
        sexes: bool
            default = False; individuals are hemraphroditic. If True,
            individuals will be male and female
        """
        ClonalHaploidWF(
            model=self.model,
            pop_size=pop_size,
            sexes=sexes,
        ).run()

    def wright_fisher_altgen(
        self,
        spo_pop_size: int,
        gam_pop_size: int,
        sexes: bool=False,  # mono/hetero terms is more consistent with others...
        ):
        """
        Generate scripts appropriate for an altered Wright-Fisher model 
        with alternation of generations. Each generation reproduces and
        experiences fitness effects.

        Parameters:
        -----------
        spo_pop_size: 
            Size of the population in number of diploids in diploid life stage
        gam_pop_size:
            Size of the population in number of haploids in the haploid life stage
        sexes: bool
            default = False; individuals are hemraphroditic. If True,
            individuals will be male and female
        """
        AltGenWF(
            model=self.model,
            spo_pop_size=spo_pop_size,
            gam_pop_size=gam_pop_size,
            sexes=sexes,
        ).run()
