initCobraToolbox
modelFileName = 'iRsp_1066_v1.5.xml';
changeCobraSolver('gurobi6', 'ALL')
dir = fileparts (which (modelFileName));
model = readCbModel(modelFileName);
ori_model = model ;

cd /home/ioannisvasilas
model.c(strcmp(model.rxns, 'BIOMASS_aerobic')) = 0.01
model = changeRxnBounds(model ,'EX_nh4_e' ,-1000, 'l');
model = changeRxnBounds(model ,'ATPPH' ,8.39, 'l');
model = changeRxnBounds(model ,'EX_o2_e' ,-1000, 'l');
model = changeRxnBounds(model ,'EX_glc__D_e' ,-1000, 'l');
model = changeRxnBounds(model ,'EX_k_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_na1_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_ca2_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_pi_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_so4_e' ,-1000 , 'l');
model = changeRxnBounds (model ,'EX_zn2_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_fe2_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_mn2_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_cobalt2_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_nac_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_thm_e' ,-1000, 'l');
model = changeRxnBounds (model ,'EX_btn_e' ,-1000, 'l');




growthRate = optimizeCbModel (model);
fprintf('The maximum growth rate is %1.2f', growthRate.f);

model = changeObjective ( model, 'DM_amd') ;
max_target = optimizeCbModel(model);
fprintf('The maximum production rate of terpene is %1.2f', max_target.f);

constrWT = struct('rxnList', {{'BIOMASS_aerobic'}}, 'rxnValues', 1, 'rxnBoundType', 'b')
constrMT = struct('rxnList', {{'BIOMASS_aerobic', 'DM_amd'}}, 'rxnValues', [0, 1], ...
                  'rxnBoundType', 'bb')


[minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(model, ...
                                                                     constrWT, constrMT);
                                                                 
disp([minFluxesW, maxFluxesW, minFluxesM, maxFluxesM]);

runID = 'TestOptForce_DM_amd';


 
disp(mustLSet)

[mustUSet, pos_mustU] = findMustU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
                                  'outputFileName', 'MustU' , 'printExcel', 1, 'printText', 1, ...
                                  'printReport', 1, 'keepInputs', 1, 'verbose', 0);



[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
    findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);

disp(mustUU);

[mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
    findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);


           
           
disp(mustLL);  
constrOpt = struct('rxnList', {{'EX_glc__D_e', 'BIOMASS_aerobi', 'DM_amd'}}, 'values', [-1000, 0, 1]');
exchangeRxns = model.rxns(cellfun(@isempty, strfind(model.rxns, 'EX_')) == 0);
excludedRxns = unique([mustUSet; mustLSet; exchangeRxns])

[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
    findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);

           
  disp(mustUU);
  
  [mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
    findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);

 disp(mustLL);
 
 [mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = ...
    findMustUL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUL', 'outputFileName', 'MustUL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);

disp(mustUL);

mustU = unique(union(mustUSet, mustUU));
mustL = unique(union(mustLSet, mustLL));
targetRxn = 'DM_amd';
biomassRxn = 'BIOMASS_aerobi';
k = 1;
nSets = 1;
constrOpt = struct('rxnList', {{'EX_glc__D_e','BIOMASS_aerobi'}}, 'values', [-1000, 0]);

[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
             'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
             'runID', runID, 'outputFolder', 'OutputsOptForce', ...
             'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
             'printReport', 1, 'keepInputs', 1, 'verbose', 1);




targetRxn = 'DM_amd';

biomassRxn = 'BIOMASS_aerobic';

k = 10;

nSets = 30;

constrOpt = struct('rxnList', {{'EX_glc__D_e' ,'BIOMASS_aerobic'}}, 'values', [-1000, 0]');



%% 
k = 2;

nSets = 30; % the number of posibility 

runID = 'TestOptForce';

excludedRxns = struct('rxnList', {{'IR05641'}}, 'typeReg','U'); 

[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
'excludedRxns', excludedRxns, ...
'runID', runID, 'outputFolder', 'OutputsOptForce', ...
'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
'printReport', 1, 'keepInputs', 1, 'verbose', 1);
%% Results
dataset(optForceSets, typeRegOptForceSets)

formula = printRxnFormula(model,model.rxns,false, true,true);

[loca,locb] = ismember(optForceSets(:,1),model.rxns);

[loca2,locb2] = ismember(optForceSets(:,2),model.rxns);

fid = fopen('OptForce_malcoa_result.txt','w+');
for i = 1: length(optForceSets)
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',optForceSets{i,1},optForceSets{i,2},typeRegOptForceSets{i,1},typeRegOptForceSets{i,2},formula{locb(i)},formula{locb2(i)});
end