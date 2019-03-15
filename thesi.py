# Find all the reactions
# Checks all the imports 
from __future__ import absolute_import
from __future__ import print_function
# # Loading the Libraries
import cobra.test


from cobra import Model, Reaction, Metabolite  # Best practise: SBML IDs
from cobra.test import test_all
from os.path import join
from cobra.flux_analysis import \
     single_gene_deletion, single_reaction_deletion, \
     double_gene_deletion, double_reaction_deletion
from time import time
from lxml import etree
import pickle
import scipy.io
import sys
from io import StringIO
import os
import os.path
import pandas
import optlang.interface
import numpy as np
import pandas as pd
import pickle 




# # # Loading the Model iRsp1140.xml from Imam et al. 2013

model = cobra.io.read_sbml_model("iRsp1066_BiGG.xml")
# Testing that the model was loaded correctly (by e.g. printing the number of reactions in the model)
print(len(model.reactions))


for r in model.reactions:
    if r.reaction:
        print (r,r.name,r.id,r.genes)

       

        

        


