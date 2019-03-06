import cobra
import numpy as np
model = cobra.io.read_sbml_model("iRsp1066_BiGG.xml")


oxygen_name = "EX_o2_e_"

glucose_name = "EX_glc_e_"

biomass_name = 'BIOMASS_aerobic' 

acetate_name = "EX_ac_e"

o2_rates = np.arange(initial_o2,final_o2 + step_size ,step_size)

glu_rates = np.arange(initial_glu,final_glu + step_size ,step_size)


values_growth = np.zeros((o2_rates.shape[0], glu_rates.shape[0]))

pfba_values_growth = np.zeros((o2_rates.shape[0], glu_rates.shape[0]))


with model as model:

    for i, o2_rate in enumerate(o2_rates):

        for j, glu_rate in enumerate(glu_rates):

            model.reactions.get_by_id(oxygen_name).bounds = (-o2_rate, -o2_rate)

            model.reactions.get_by_id(glucose_name).bounds = (-glu_rate, -glu_rate)

            

try:

                fba_sol = model.optimize()

                values_growth[i, j] = fba_sol.objective_value

except:

                values_growth[i, j] = 0


try:

                solution = cobra.flux_analysis.pfba(model)

                pfba_values_growth[i, j] = solution.fluxes[biomass_name]

except:

                pfba_values_growth[i, j] = 0


fig, ax = plt.subplots(figsize=(10,10))

pos = ax.imshow(ijn1411_values_growth, interpolation='none')

cbar = fig.colorbar(pos, ax=ax, fraction=0.034, pad=0.04)

cbar.set_label('Biomass production rate [mmol per gDW per hr]', rotation=270, labelpad=20, fontsize = 14)


fig.suptitle('Phenotypic Phase Plane of ', fontsize = 20)

ax.set_xlabel('Glucose Uptake Rate [mmol per gDW per hr]', fontsize = 14)

ax.set_ylabel('Oxygen Uptake Rate [mmol per gDW per hr]', fontsize = 14)


ax.imshow(ijn1411_values_growth, interpolation='none', extent=[initial_glu, final_glu,initial_o2, final_o2], origin='lower')

ax.set_aspect(0.2) # you may also use am.imshow(..., aspect="auto") to restore the aspect ratio

fig.savefig('Phenotypic Phase .png', dpi=100) 

                



