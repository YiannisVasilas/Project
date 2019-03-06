import cobra
from cobra.flux_analysis import flux_variability_analysis
from copy import deepcopy
#from pulp import LpProblem, LpMaximize, LpMinimize, LpVariable, LpAffineExpression, \
                #solvers, LpContinuous, LpBinary, LpStatusOptimal, lpSum, LpStatus

import pandas as pd
#before;
#from cobra.core.Solution import Solution
from cobra.core.solution import Solution  #I replaced by solution in small letters

#TO DO
# Try accessing Cplex directly as alternative to OptSlope code 
# Replace 'yes'/'no' tilt by True/False       
# KO_cand is sometimes '' and sometimes None. Standardize.
# Add integer cuts
# Allow single reaction (non-list) to be added in OptAndRob

def OptAndRob(model,targets,numDel,KO_cand,verbose=False): # add again the KO_cand=''
    ### PURPOSE
    #   Performs OptKnock and RobustKnock for a model for any number of targets
    ### EXAMPLE CALL
    #   Opt_R,Rob_R = OptAndRob(model,ExchangeReactions,1)
    ### INPUT
    #   model:    CobraPy model structure
    #   targets:  A list of CobraPy reaction IDs that are to be optimised (individually)
    #   numDel:    The number of KOs for OptKnock and RobustKnock
    ### OPTIONAL INPUT
    #   KO_cand:  A list of reaction IDs for reactions that can be knocked out. Default: all reactions in model
    #     verbose:  A boolean indicator of whether intermediate results should be printed to screen. Default: False.
    ### OUTPUT
    #   Opt_R:    OptKnock results
    #   Rob_R:    RobustKnock results
    
    #If no candidate target reactions are given set all reactions in the model as targets
    if KO_cand=='':
        KO_cand=[r.id for r in model.reactions]

    #Create output structures
    Opt_R = {}
    Rob_R = {}

    #Determine number of targets
    if type(targets)==str:
        targets = [targets]
    totalNum = len(targets) #If single reaction added this does not work!!
    
    for i,target in enumerate(targets):
        if verbose:
            print("At number "+str(i+1)+" of "+str(totalNum)+" targets")

        #OptKnock
        okmodel = Run_OptKnock(model,target,numDel,verbose=False,integralityTol=0.,KO_cand=KO_cand)
        if hasattr(okmodel.solution,'KOs'): #I replace okmodel.solution  by optimize()
            KOs = [r.id for r in okmodel.solution.KOs]
            minval,maxval,gr = evaluateKO(model,target,KOs) #Use FBA to double-check solution
        else:
            KOs = []
            minval = 0
            maxval = 0
            gr = 0
        Opt_R[target] = {'KOs':KOs,'minimum':minval,'maximum':maxval,'growthRate':gr}

        #RobustKnock
        okmodel = Run_OptKnock(model,target,numDel,tilt=True,verbose=False,integralityTol=0.,KO_cand=KO_cand)
        if hasattr(okmodel.solution,'KOs'):
            KOs = [r.id for r in okmodel.solution.KOs]
            minval,maxval,gr = evaluateKO(model,target,KOs) #Use FBA to double-check solution
        else:
            KOs = []
            minval = 0
            maxval = 0
            gr = 0
        Rob_R[target] = {'KOs':KOs,'minimum':minval,'maximum':maxval,'growthRate':gr}
        
    return Opt_R,Rob_R

def Run_OptKnock(model,targetRxn,numDel,objRxn='',objFraction=0.2,verbose=True,tilt=False,tiltValue=0.00001,epgap=0.0001,integralityTol=0.,KO_cand = None,numSol=1):
    ### PURPOSE
    #   Runs the OptKnock algorithm on the model.
    ### EXAMPLE CALL
    #   results = Run_OptKnock(model,'EX_Succinate',3)
    ### INPUT
    #   model:            CobraPy model structure
    #    targetRxn:        Reaction ID of the reaction that is to be optimised.
    #     numDel:            The maximum number of allowed KOs.
    ### OPTIONAL INPUT
    #   objRxn:         Reaction ID corresponding to the cellular objective. By default OptKnock will check the cobra model for the objective reaction.
    #    objFraction:    The minimal fraction of normal growth that is maintained. 
    #    verbose:         A boolean indicator of whether intermediate results should be printed to screen. Default: False.
    #     tilt:            Whether or not to tilt the objective. 'no' corresponds to the original optimistic OptKnock, 'yes' corresponds to a pessimistic version more similar to RobustKnock.
    #    tiltValue:        The size of the (negative!) tilt of the targetRxn to the objective function.
    #     epgap:            Acceptable difference between the final result and the dual result. 
    #     integralityTol: Maximal deviation of a integer variable from an integer.
    #     KO_cand:        The reaction IDs of reactions that may be knocked out. Default: all reactions
    #     numSol:         Number of alternative solutions that are desired. These are obtained using Integer cuts.
    ### OUTPUT
    #   okmodel:         OptKnock model structure with solution.

    #Process input
    assert targetRxn in [r.id for r in model.reactions], "targetRxn is not in the model: %r" % targetRxn
    assert type(numDel) is int, "numDel is not an integer: %r" % numDel
    assert objRxn in [r.id for r in model.reactions] or not objRxn, "targetRxn is not in the model: %r" % objRxn
    assert 0<=objFraction<=1, "The objective fraction has to be between 0 and 1. The value is: %r" % objFraction
    assert verbose==True or verbose==False, "Verbose has to be True or False." % objFraction
    assert tilt==False or tilt==True,"Tilt has to be either yes or no."
    assert epgap>=0,"epgap has to be larger or equal to 0."
    assert integralityTol>=0,"integralityTol has to be larger or equal to 0."

    #Set default objective reaction
    if not objRxn:
        objRxn = [rxn.id for rxn in model.reactions if rxn.objective_coefficient>0]
        assert len(objRxn)==1,'There should be exactly one objective reaction. The number of objective reactions is: %d' % len(objRxn) #Code does not work with multiple -> terminate.
        objRxn = objRxn[0]

    #Initialize model
    okmodel = OptKnock(model,verbose) 

    #Set minimal growth rate
    if objFraction > 0:
        okmodel.model.optimize()
        okmodel.model.reactions.get_by_id(objRxn).lower_bound = okmodel.model.optimize().objective_value*objFraction # i changed model.solution.f by model.Solution.objective_value

    #Identify KO candidates
    if KO_cand is not None:
        KO_cand = [r for r in okmodel.model.reactions if r.id in KO_cand]

    #Prepare problem
    okmodel.prepare_optknock(targetRxn,ko_candidates=KO_cand,num_deletions=numDel,tilt=tilt,tiltValue=tiltValue) #I introduce self

    #Set parameters in solver
    # epgap
    okmodel.prob.solver.epgap = epgap
    # integrality tolerance
    okmodel.prob.solver.integralityTol = integralityTol

    #Solve problem
    okmodel.solve()

    #Print results, if verbose
    if verbose:
        PrintResults(okmodel,targetRxn)

    return okmodel

#This function serves to evaluate the proposed KOs in terms of minimal and maximal production values
def evaluateKO(model,targetRxn,KO_IDs):
    ### PURPOSE
    #   Predicts the minimal and maximal fluxes of the target reaction and of the biomass reaction based on KOs
    ### EXAMPLE CALL
    #   evaluateKO(model,results)
    ### INPUT
    #   model:      CobraPy model structure
    #     targetRxn:    Reaction ID of the reaction that is to be optimised.
    #     KO_IDs:        Reaction IDs of the reactions that have been knocked out.
    ### OPTIONAL INPUT
    ### OUTPUT

    #Create temporary model - make no changes to original model
    tmpmodel = deepcopy(model)

    #Mutant
    tmpmodel = block_reactions(tmpmodel,KO_IDs,cumulative=True)
    tmpmodel.optimize()
    growthrate = tmpmodel.optimize().objective_value # I change tmpmodel.solution.f by objective_value()

    if growthrate > 10**-6: #KOs are not lethal
        #Minimal and maximal production during maximum growth
        #Before:
        #FVA = cb.flux_analysis.variability.flux_variability_analysis(tmpmodel,fraction_of_optimum=1,reaction_list=[tmpmodel.reactions.get_by_id(targetRxn)])
        FVA = cobra.flux_analysis.flux_variability_analysis(tmpmodel,reaction_list=[tmpmodel.reactions.get_by_id(targetRxn)],fraction_of_optimum=1) #I modify the cb and put cobra. Remove the .variability		
        minval=FVA.minimum[targetRxn]
        maxval=FVA.maximum[targetRxn]
		#Before it was like this:
        # minval = FVA[targetRxn]['minimum']
        # maxval = FVA[targetRxn]['maximum']
    else: #Proposed KOs are lethal, correct solution was not found.
        minval = float('nan')
        maxval = float('nan')
        growthrate = 0
    
    return minval,maxval,growthrate

def PrintResults(okmodel,targetRxn):
    ### PURPOSE
    #   Prints the (tilted) Optknock results as well as the results of an independent evaluation 
    #   of the KOs using FBA. This is required as OptKnock solutions are occasionaly incorrect due to numerical issues
    ### EXAMPLE CALL
    #   PrintResults(okmodel)
    ### INPUT
    #   okmodel:     An Optknock model object.
    #     targetRxn:    Reaction ID of the reaction that is to be optimised.
    ### OPTIONAL INPUT
    ### OUTPUT
    #     Prints the (tilted) OptKnock results to screen.

    print('\nObjective: %6.3f' %okmodel.prob.objective.value())
    print('Biomass rate %6.3f:'%okmodel.optimize().x_dict[okmodel.r_biomass])
    print("Sum of mu : %6.3f" % sum(okmodel.optimize().mu))

    print '\nKnockouts:flux'
    for r in okmodel.solution.KOs:
        print('%s:%6.10f' % (r.id, okmodel.solution.x_dict[r]))

    print '\nNon-binary y''s'
    for r in okmodel.model.reactions:
        val = okmodel.solution.y_dict[r]
        if val!=0 and val!=1:
            print('%s:%6.10f' % (r.id, okmodel.solution.x_dict[r]))

    #Use FBA to double-check solution
    KOs = [r.id for r in okmodel.solution.KOs]
    minval,maxval,gr = evaluateKO(okmodel.model,targetRxn,KOs) 
    print '\nResults of independent FBA evaluation'
    print('Biomass rate %6.3f:'%gr)
    print('Minimal target rate %6.3f:'%minval)
    print('Maximal target flux %6.3f:'%maxval)

