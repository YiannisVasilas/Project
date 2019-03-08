# coding=utf-8 # because of copy-pasting non-ASCII conversion from web
# keep this block always ON
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





# # # Loading the Model iRsp1140.xml from Imam et al. 2013

model = cobra.io.read_sbml_model("iRsp1066_BiGG.xml")
# Testing that the model was loaded correctly (by e.g. printing the number of reactions in the model)
print(len(model.reactions))



# # # In Silico Growth - Experiment 1
# # Step 1: values in vmin and vmax were set to −100 and 100 mmol/g DW h for reversible reactions
# 1.1 Find reversible reactions
# 1.2 Change reversible reactions allowable fluxes to -100 to 100
# # Step 2: Figure out which are the media components essential for cell growth
# 2.1 change boundary condition of compounds where there are measured rates available to those rates
# 2.2 change boundary conditions to -100 and 100
# # Step 3: Adjust non-growth associated ATP Maintenance (NGAM) to 8.39 mmol/g DW h
# ??? I guess I add a reaction which utilizes ATP -> ADP + Pi and I set the boundary conditions to: 0 (lower) 8.39 (upper)?!
# # Step 4: Adjust objective function
# change to what? Imam2013 does not give an indication what it should be changed to. So I assumed aerobic growth at this moment (27th Dec)


# # Step 1:
# adjusting the lower and upper boundaries to -100 and 100
for r in model.reactions:
	if r.reversibility: # step 1.1
		r.lower_bound = -100 # step 1.2
		r.upper_bound = 100 # step 1.2
# # checking if lower and upper boundaries were adjusted to -100 and 100
for r in model.reactions:	
 	if r.reversibility:
 		#print(r.id,r.name, r.lower_bound, r.upper_bound)


# # Step 2:
# Step 2.1:


# ???
# Step 2.2:


print (lower) print (upper)
# # Step 3:
# ??? I guess I add a reaction which utilizes ATP -> ADP + Pi and I set the boundary conditions to: 0 (lower) 8.39 (upper)?!


# # Step 4:
# # Changing objective function to aerobic biomass optimization

biomass_rxn= model.reactions.get_by_any('BIOMASS_aerobic')
from cobra.util.solver import linear_reaction_coefficients
                                       
model.objective = ("BIOMASS_aerobic")
model.reactions.get_by_any("BIOMASS_aerobic").upper_bound=-1000
linear_reaction_coefficients(model)



                              # optimize for aerobic
print (model.optimize())
print (model.summary())


# # Checking all reactions
for r in model.reactions:
 	print(r, r.name, r.lower_bound, r.upper_bound)









# # boundary conditions











print (cobra.manipulation.modify.initialize_growth_medium(model, the_medium='MgM', external_boundary_compartment='e', external_boundary_reactions=None, reaction_lower_bound=0.0, reaction_upper_bound=1000.0, irreversible=False, reactions_to_disable=None))

boundid = []
for r in model.reactions:
 	if (r.boundary and not(r.reversibility)):
 		boundid.append(r.id)
print (len(boundid))

print (r.reversibility)
print (not(r.reversibility))
print (not(1 == 1))

boundreact = []
filter ((x for x in boundid if x in 
boundid in model.reactions.boundreact))


print(model.optimize())






# # test which metabolites are already in the model and which ones are not 


for r in model.reactions:
        if r.boundary:
                print(r, r.name, r.lower_bound, r.upper_bound)
                r.lower_bound = -100
print(model.optimize())
print(model.summary())



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




