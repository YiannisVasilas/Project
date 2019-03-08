import cobra
import os

model = cobra.io.read_sbml_model('iRsp1140rc.xml')
def flux_uptake(model):
    model = cobra.io.load_sbml_model()
    if model:
                    model = model.reactions.get_by_id(uptake_id)
                    model_reac.lower_bound = -100.0
                    model_reac.upper_bound = -100.0
                    print ("Uptake reaction {} is fixed to {}.".format(uptake_id, uptake_reac.upper_bound))
                    print (flux_uptake)

            
    
