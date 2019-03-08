
# coding: utf-8

# # iRsp1066 OptGene 2

# ## Loading the Libraries

# In[3]:

from __future__ import division # for being able to divide by "/" character

## came
import cameo
from cameo import fba # for doing a flux balance analysis
#from cameo.exceptions import Infeasible
from cameo import load_model
from cameo import models
## strain design
from cameo.strain_design.deterministic.linear_programming import OptKnock
from cameo.strain_design.heuristic.evolutionary_based import OptGene
## for phaseplane
from cameo.visualization import plotting
from cameo import phenotypic_phase_plane

# for FVA
from cameo import flux_variability_analysis

# for phaseplane
from cameo.visualization import plotting
from cameo import phenotypic_phase_plane

import cobra
from cobra import Reaction ,Model, Metabolite

import numpy
import pandas

## Debugging
from pprint import pprint # for showing details of the model with pprint(vars(model))

import string # to replace a string with another string


# ## Loading the Model iRsp1066.xml

# In[4]:

# model = cobra.io.read_sbml_model("/home/mario/Experiments/StrainConstruction/iRsp265_updated_MVA.xml") #cobra
model = load_model("iRsp1066.xml")
BOF = "RXN1391"

NrKnockOuts = 2

# RUN = "OK%s" % str(NrKnockOuts)
RUN = "OG%s" % str(NrKnockOuts)      # no growth coupling
# RUN = "OG%s_gc" % str(NrKnockOuts)   # growth coupled run


# In[ ]:

print( RUN)


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
    PO4         = copymodel.reactions.get_by_id("RXN0213")        # Phosphate (PO4)
    SO4         = copymodel.reactions.get_by_id("RXN0196")        # Sulfate (SO4)
    Mg2         = copymodel.reactions.get_by_id("RXN1158")        # Magnesium (Mg2+)
    Fe2         = copymodel.reactions.get_by_id("RXN0192")        # Iron (Fe2+)
    Nicotinate  = copymodel.reactions.get_by_id("RXN1329")        # Nicotinate
    Thiamine    = copymodel.reactions.get_by_id("RXN1354")        # Thiamine
    Biotin      = copymodel.reactions.get_by_id("RXN1326")        # Biotin

    MediaComponents = [PO4, SO4, Fe2, Nicotinate, Thiamine, Biotin, Mg2]

    for r in MediaComponents:
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
    Glucose.lower_bound = -4
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
model.reactions.get_by_id("RXN1922").lower_bound = 0.0
model.reactions.get_by_id("RXN1922").upper_bound = 1000.0

## setting objective function
# model.change_objective("RXN1922") # Amorphadiene
# Do not set Amorphadiene as the objective. The "system" objective should remain biomass.

## setting minimum BOF lower growth bound
# using ~30% of the max growth rate under wt conditions 0.09
model.reactions.get_by_id(BOF).lower_bound = 0.03
model.reactions.get_by_id(BOF).upper_bound = 1000.0


# ## OptGene

# In[8]:

optgene = OptGene(model)

result_OG = optgene.run(target = "RXN1922", 
                     # target = "EX_ac_e",
                     biomass = BOF,
                     substrate = "CPD1074",
                     max_knockouts = NrKnockOuts,
                     variable_size = True,
                     simulation_method = "fba",
                     growth_coupled = False,
                     max_evaluations = 20000,
                     population_size = 200,
                     time_machine = None,
                     max_results = 50,
                     use_nullspace_simplification = True,
                        plot = False)
print (model.summary())



OptGeneResult = pandas.DataFrame(result_OG.data_frame)
OptGeneResult.to_csv(path_or_buf = "/home/Results/%s.csv" % RUN, sep = ";", decimal = ",")

