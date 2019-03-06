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

from cameo.strain_design.deterministic.linear_programming import OptKnock
from cameo.strain_design.heuristic.evolutionary_based import OptGene
## for phaseplane
from cameo.visualization import plotting
from cameo import phenotypic_phase_plane

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
import gurobipy5
model = cobra.io.read_sbml_model("iRsp1066_BiGG.xml")

print(len(model.reactions))

#optgene = OptGene(model)
#result = optgene.run(target=model.reactions.DM_amd,
                     #biomass=model.reactions.BIOMASS_aerobic,
                     #substrate=model.metabolites.glc__D_e,
                     #max_evaluations=50000,
                     #plot=False)
                     
#print (result)


optknock = OptKnock(model, fraction_of_optimum=0.1)
result = optknock.run(max_knockouts=1, target= "DM_ipdp", biomass="BIOMASS_aerobic")
print (result)
