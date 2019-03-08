# coding=utf-8 # because of copy-pasting non-ASCII conversion from web
# keep this block always ON

# # Loading the Libraries
import cobra.test
#from __future__ import print_function
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
#import cStringIO
import os
import os.path
import pandas


# Loading the Model iRsp1140.xml from Imam et al. 2013
model = cobra.io.read_sbml_model("Rsp1140.xml")
# Testing that the model was loaded correctly (by e.g. printing the number of reactions in the model)
print(len(model.reactions))


# # # In Silico Growth - Experiment 1
# # Step 1: values in vmin and vmax were set to −100 and 100 mmol/g DW h for reversible reactions
# Find reversible reactions
# Change reversible reactions allowable fluxes to -100 to 100



print (cobra.flux_analysis.essentiality.assess_medium_component_essentiality(model, the_components=None, the_medium=None, medium_compartment='e', solver='glpk', the_condition=None, method='fba'))











print(model.optimize())



model.change_objective("RXN1391")
print(model.optimize())
#test which metabolites are already in the model and which ones are not 


for r in model.reactions:
 	if r.boundary:
 		print(r, r.name, r.lower_bound, r.upper_bound)
 		lower_bound = -100
 		print(model.optimize())





'''Documentation with Maria Talk'''
# which carbon and nitrogen sources are allowing growth? from the paper
## Model
# evaluate in silico carbon and nitrogen sources (from the paper)
# add reactions of the mevalonate pathway (email from Frank)

# are there genes in the model which are not in the industrial strain?
# if they are isozymes, then deleting would not change anything?

## Annotation
# how many genes do they have in comon? -> venn diagramm

## eLabNotebook

# write more emails/ ask if I am stuck
# send a short email of how I am doing with the model after the exam





