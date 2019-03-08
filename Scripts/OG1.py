
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

## Debugging
from pprint import pprint # for showing details of the model with pprint(vars(model))

import string # to replace a string with another string


# ## Loading the Model iRsp1066.xml

# In[4]:

model = cobra.io.read_sbml_model("iRsp1066.xml") #cobra
#model = load_model("iRsp1066.xml")
BOF = "RXN1391"







# ## Initiate Growth Media

# In[5]:

def initiateMedia(copymodel = model):

    ## Global Variables
    LB = -1000
    UB = 1000


    ## Adjusting values of v<sub>min</sub> and v<sub>max</sub>
    for r in copymodel.reactions:
        if r.reversibility: # Finding reversible reactions
            r.lower_bound = LB
            r.upper_bound = UB


    ## Adjusting Boundary Conditions
    for r in copymodel.reactions:
        if r.boundary:
            r.lower_bound = 0
            r.upper_bound = UB

    for r in copymodel.reactions:
        if (not r.boundary and not r.reversibility):
            r.lower_bound = 0
            r.upper_bound = UB


    ## Adjusting all Boundary Conditions to ± LB/ UB
    for r in copymodel.reactions:
        if r.lower_bound == -1000:
            r.lower_bound = LB
        if r.upper_bound == 1000:
            r.upper_bound = UB


    ## Adjusting the Non-Growth Associated ATP Maintenance Reaction
    ATPM = copymodel.reactions.get_by_id("RXN0765")
    ATPM.lower_bound = 8.39
    ATPM.upper_bound = 8.39


    ## Adjusting the Media Components
    ## Adjusting the Uptake Rates of Media Components
    Potassium = model.reactions.get_by_id('RXN1167')
    Aspartic_Acid = model.reactions.get_by_id('RXN1474')
    Sodium = model.reactions.get_by_id('RXN0220')
    Calcium = model.reactions.get_by_id('RXN0221')
    Phospate = model.reactions.get_by_id('RXN0213')
    Sulfate = model.reactions.get_by_id('RXN0196')
    Zinc = model.reactions.get_by_id('RXN0190')
    Iron = model.reactions.get_by_id('RXN0192')
    Magnesium = model.reactions.get_by_id ('RXN1158')
    Copper = model.reactions.get_by_id ('RXN0188')
    #Cobalt = model.reactions.get_by_id ('EX_cobalt2_e')
    Nicotin = model.reactions.get_by_id('RXN1329')
    Thiamine =model.reactions.get_by_id('RXN1354')
    Biotin = model.reactions.get_by_id('RXN1326')

    Media_Components= [Potassium, Aspartic_Acid,Phospate, Sodium, Calcium ,Zinc,Iron, Magnesium, Copper, Nicotin ,Thiamine, Biotin]



    

    for r in Media_Components:
        r.lower_bound = LB
        r.upper_bound = UB


    ## Aerobic Growth Conditions
    O2          = copymodel.reactions.get_by_id("RXN0223")        # Oxygen (O2)
    # H2O         = copymodel.reactions.get_by_id("RXN0217")        # H2O
    # H             = copymodel.reactions.get_by_id("RXN1395")      # H+
    
    SmallMolecules = [O2] #, H2O]

    for r in SmallMolecules:
        r.lower_bound = LB
        r.upper_bound = UB
        
    ## H+
    # H             = copymodel.reactions.get_by_id("RXN1395")      # H+
    # H.lower_bound = -2.0    
    # H.upper_bound = 1000.0
    # Petersen et al. (2012) --> doi: 10.1073/pnas.1202799109
    # ratio of H+/ATP = 4+-0.3
    # shadow_price of ATP (CDP0002) was tested under varying H2O and H+ levels
    # shadow price of ATP was always: -0.479 --> -0.5*4 = -2
    # thus flux of H+ was constrained to +-2


    ## Adjusting Nitrogen Uptake Rates
    NH4 = copymodel.reactions.get_by_id("RXN0219")
    NH4.lower_bound = -1 # (should be -1. But -10 makes sure that the carbon source is the limiting factor)
    NH4.upper_bound = 0


    ## Carbon Source Uptake Rate
    Glucose = copymodel.reactions.get_by_id("RXN1101") # glucose: RXN1101/ glycerol: RXN0209
    Glucose.lower_bound = -100
    Glucose.upper_bound = 0

    ## Switching off all Demand reactions
    demandRXNs = numpy.array([[rxn.id, rxn.lower_bound, rxn.upper_bound, rxn.name] for rxn in model.reactions if "demand" in str(rxn.name).lower()])

    for rxn in model.reactions:
        if rxn.id in demandRXNs[:,0]:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

    ## Change model objective
    copymodel.objective=("RXN1391") # aerobic
    
    return copymodel;


# ## Example
# copymodel = model.copy()
# copymodel = initiateMedia(copymodel)
# copymodel.solve()


# In[6]:

model = initiateMedia(model)


# In[7]:

## switching on Amorphadiene Demand
model.reactions.get_by_id("RXN0668").lower_bound = 0
model.reactions.get_by_id("RXN0668").upper_bound = 1000


## setting objective function
# model.change_objective("RXN1922") # Amorphadiene
# Do not set Amorphadiene as the objective. The "system" objective should remain biomass.

## setting minimum BOF lower growth bound
# using ~30% of the max growth rate under wt conditions 0.09
model.reactions.get_by_id(BOF).lower_bound = 0.03
model.reactions.get_by_id(BOF).upper_bound = 1000.0
#print (model.slim_optimize())
print(model.medium)
#from cobra .medium import minimal_medium
#max_growth= model.slim_optimize()
#print (minimal_medium(model,max_growth))
#print (len (minimal_medium(model,max_growth)))
cobra.io.write_sbml_model( model, "test_cobra.xml", use_fbc_package=False)
#en=production_envelope  (model,['EX_o2_e'],objective= 'BIOMASS_aerobic', carbon_sources= 'EX_glc__D_e',points = 100)
en =np.array(production_envelope  (model,['RXN0223'],objective= 'RXN1391', carbon_sources= 'RXN1101',points = 100))
np.savetxt('En2.csv',en,fmt='%s')
print (en)



