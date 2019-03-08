
from __future__ import division #
from __future__ import absolute_import
import cobra.test
import random, string, time

from cobra import Model, Reaction, Metabolite  
from cobra.test import test_all
from os.path import join
from cobra.flux_analysis import \
     single_gene_deletion, single_reaction_deletion, \
     double_gene_deletion, double_reaction_deletion, production_envelope
from time import time
from lxml import etree


## came


## strain design
from cameo.strain_design.deterministic.linear_programming import OptKnock
from cameo.strain_design.heuristic.evolutionary_based import OptGene

import cobra
from cobra import Reaction , Model, Metabolite

import numpy
import pandas

## Debugging
from pprint import pprint # for showing details of the model with pprint(vars(model))

import string # to replace a string with another string

import pickle
import scipy.io
import sys
from io import StringIO
import os
import os.path
import pandas
import numpy as np
from cobra.exceptions import OptimizationError, SolverNotFound,\
    OPTLANG_TO_EXCEPTIONS_DICT,Infeasible
from cobra.util.context import get_context
import optlang.interface






### Load mdel and check if loaded 

model = cobra.io.read_sbml_model("iRsp1066_BiGG.xml")

print(len(model.reactions))

biomass = 'BIOMASS_aerobic' 




        # Carbon and Nitrogen list for the revaluation
'''Carbon_Exchange_List={'EX_glc(e)','EX_pyr(e)','EX_glyc(e)','EX_ac(e)','EX_rib_D(e)' ,'EX_fru(e)','EX_sbt_D(e)','EX_gal(e)','EX_lcts(e)','EX_man(e)'}
Nitrogen_Exchange_List={'EX_nh4(e)','EX_ade(e)','EX_cytd(e)','EX_ptrc(e)','EX_gly(e)','EX_ala_L(e)','EX_gln_L(e)','EX_arg_L(e)','EX_gam(e)','EX_4abut(e)'}'''


# Adjusting all the reversible reactions
for r in model.reactions:
        if r.reversibility:
                r.lower_bound = -1000
                r.upper_bound = 1000

#for r in model.reactions:       
        #if r.reversibility:
                
                #print(r.id,r.name, r.lower_bound, r.upper_bound)

# Set boundary conditions
for r in model.reactions:
        if r.boundary:
                r.lower_bound = 0
                r.upper_bound = 1000

for r in model.reactions:
        if (not r.boundary and not r.reversibility):
                r.lower_bound = 0
                r.upper_bound = 1000
                                         



'''for r in model.reactions:
        if r.lower_bound == -1000:
                r.lower_bound = -1000
                if r.upper_bound == 1000:
                        r.upper_bound = 1000'''
                

        
# Non Grotwth  ATP Maintenance Reaction
ATPM = model.reactions.get_by_id('ATPPH')

ATPM.lower_bound = 8.39
ATPM.upper_bound = 8.39



#Media component

Potassium = model.reactions.get_by_id('EX_k_e' )
Aspartic_Acid = model.reactions.get_by_id('EX_gly_asp__L_e')
Sodium = model.reactions.get_by_id( 'EX_na1_e')
Calcium = model.reactions.get_by_id( 'EX_ca2_e')
Phosphate = model.reactions.get_by_id('EX_pi_e')
Sulfate = model.reactions.get_by_id('EX_so4_e')
Zinc = model.reactions.get_by_id('EX_zn2_e')
Iron = model.reactions.get_by_id( 'EX_fe2_e')
Magnesium = model.reactions.get_by_id ('EX_mn2_e')
Copper = model.reactions.get_by_id ('EX_cu2_e')
Cobalt = model.reactions.get_by_id ('EX_cobalt2_e')
Nicotin = model.reactions.get_by_id('EX_nac_e')
Thiamine =model.reactions.get_by_id('EX_thm_e')
Biotin = model.reactions.get_by_id('EX_btn_e')

Media_Componets= [Potassium,Aspartic_Acid,Sodium,Calcium,Phosphate,Sulfate,Zinc,Iron,Magnesium,Copper,Cobalt,Nicotin ,Thiamine, Biotin]

for r in Media_Componets:
        r.lower_bound = -10
        r.upper_bound = 1000


#aerobic growth and nitrogen uptake
#H2 = model.reactions.get_by_id("EX_h2_e")
#H2O = model.reactions.get_by_id('EX_h2o_e')
#H2.lower_bound = -1000
#H2.upper_bound = 1000
#H2O.lower_bound = -1000
#H2O.upper_bound = 1000

O2 = model.reactions.get_by_id("EX_o2_e")
O2.lower_bound = -1000
O2.upper_bound = 1000

NH4 = model.reactions.get_by_id("EX_nh4_e")
NH4.lower_bound = -1
NH4.upper_bound = 0





# Carbon source

Glucose = model.reactions.get_by_id('EX_glc__D_e')
Glucose.lower_bound = -4
Glucose.upper_bound = 0

# demand reactions off
demandRXNs = numpy.array([[rxn.id, rxn.lower_bound, rxn.upper_bound, rxn.name] for rxn in model.reactions if "demand" in str(rxn.name).lower()])

for rxn in model.reactions:
        if rxn.id in demandRXNs[:,0]:
                rxn.lower_bound = 0
                rxn.upper_bound = 0

'''Carbon_source = model.exchanges
for r in Carbon_source:
        r.lower_bound = -100
        r.upper_bound = 14'''



#objective
                                       
model.objective = ("BIOMASS_aerobic")                            
model.reactions.get_by_id("DM_amd").lower_bound = 0
model.reactions.get_by_id("DM_amd").upper_bound = 1000
model.reactions.get_by_id('BIOMASS_aerobic').lower_bound = 0.2
model.reactions.get_by_id('BIOMASS_aerobic').upper_bound = 1000
model.reactions.get_by_id('EDA').lower_bound = 0
model.reactions.get_by_id('EDA').upper_bound = 0
print (model.optimize())
print (model.summary())

#print (model.slim_optimize())

#from cobra .medium import minimal_medium
#max_growth= model.slim_optimize()
#print (minimal_medium(model,max_growth))
#print (len (minimal_medium(model,max_growth)))
cobra.io.write_sbml_model( model, "EDA.xml", use_fbc_package=False)
#en=production_envelope  (model,['EX_o2_e'],objective= 'BIOMASS_aerobic', carbon_sources= 'EX_glc__D_e',points = 1)
en =np.array(production_envelope  (model,['DM_amd',"EX_o2_e"],objective= 'BIOMASS_aerobic', carbon_sources= 'EX_glc__D_e',points = 1000))
np.savetxt('EDA.csv',en,fmt='%s')
print (en)


#OptKnock


'''KO= OptKnock (model)

result_OK = optknock.run(max_knockouts = NrKnockOuts,
                      target        = "RXN0569",
                      biomass       = biomass,
                      max_results   = 50)

OptKnockResult = pandas.DataFrame(result_OK.data_frame)
OptKnockResult.to_csv(path_or_buf = "/home//%s.csv" % RUN, sep = ";", decimal = ",")'''

































