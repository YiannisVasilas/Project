import cobra as cb
import sys
import time
import shutil
from pathlib2 import Path
import pandas
pandas.options.display.max_rows=350
import libsbml
import escher
import json
from escher import Builder
from cobra import Model, Reaction, Metabolite
from cobra import flux_analysis
import random



model = cb.io.read_sbml_model('modelo_cur_wt_r27.xml')

database = cb.io.read_sbml_model("REFBiGG.xml")

randomPool = []
for i in range(len(database.reactions)):
    randomPool.append(i)

#reference_db , to add the reactions when applying randomization
reference_db = cb.io.read_sbml_model("EmptyModel.xml")

model.reactions.EX_glc.lower_bound=-1
model.reactions.EX_glc.lower_bound=-1
while (len(randomPool) > 0):
    rand = random.randrange(len(randomPool))
    randd = randomPool[random.randrange(len(randomPool))]
    reference_db.add_reaction(database.reactions[randd])
    randomPool.remove(randd)

#Filter reference database
Lib_1.delete_overlap(model, reference_db) 
#CO2_g.Translate_to_MNXRef(reference_db,'bigg')
model.objective='EX_curcumin'
modd = model.optimize() 
print modd.objective_value

#Combine model and reference database and initiate MILP
gf = Lib_1.GapFiller(model,universal=reference_db, penalties = dict(universal=1, exchange=100, demand=1))
#gf = GapFiller(model,universal=reference_db)


modd = model.optimize()
print modd.objective_value

for i in range(1,4): #it was written (6,7)
    gf_out = gf.fill(iterations=50, num_reactions=i)
    for sol in gf_out:
        ### Write results to file
        #Open output file
        output = open("output_gap_cur_r27.txt","a") 
    
        #Determine Identifier
        output.write('ID_'+str(New_ID_counter)+'_r'+str(i)+'\t') #ID
        New_ID_counter += 1
        
        #output.write("Ecoli_core_aerobic"+'\t') #Model_ID
        #output.write("ref_known_pathways"+'\t') #ReferenceDB_ID
        output.write(','.join([r.id for r in sol['reactions']])+'\t') #Reactions
        output.write(str(sol['growth_rate'])+'\n') #Growth rate
        #output.write(Aerobic+'\n') #Aerobic
        
        output.close()
