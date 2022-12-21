Code of the Regulated Ecosystem Model (REcoM) for coupling with the Finite Element Sea Ice-Ocean Model 1.4 (FESOM1.4)
------------------------------------------------------------------------------------------------------------------------

Code descriptions in Hauck et al. 2013 (doi:10.1002/2013GB004600), Karakus et al. 2021 (doi:10.1029/2021JC017315), and Seifert et al. 2022 (doi:10.1525/elementa.2021.00104).

The code was used in the simulations of Seifert et al. 2022 (doi:10.1525/elementa.2021.00104) and referred to in OceanNETs deliverable 4.5 (*).



Main files that were adapted to accomodate coccolithophores and CO2 dependencies (see description of code changes in Seifert et al. 2022):

Coccolithophores, CO2 dependencies, calcification, calcite dissolution:
recom.f90

Initialisation of coccolithophores:
recom_init.f90

Adaptation of mocsy routine (Orr & Epitalon 2015, doi:10.5194/gmd-8-485-2015) for computing the 3D carbonate system and CO2 dependencies:
- buffesm.f90
- gasx.f90
- vars.f90

Namelist to define parameters:
namelist.recom
-> flags are still set as they have been used in the latest simulation. See Seifert et al. 2022 for more details on the modes (true/false) of the flags for the individual simulations.



All other files are necessary to run REcoM coupled to FESOM1.4 (Schourup-Kristensen et al. 2014, doi:10.5194/gmd-7-2769-2014; Wang et al. 2014, doi:10.5194/gmd-7-663-2014)




For questions, please contact Miriam Seifert (miriam.seifert@awi.de) or Judith Hauck (judith.hauck@awi.de).


Bremerhaven (Germany), 19 December 2022


(*) European Union’s Horizon 2020 research and innovation program under grant agree- ment number 869357 (project OceanNETs: Ocean-based Negative Emission Technologies—analyzing the feasibility, risks, and co-benefits of ocean-based negative emission technologies for stabilizing the climate). The work reflects only the authors' views; the European Commission and their executive agency are not responsible for any use that may be made of the information the work contains.
