clear all
clc
initCobraToolbox

tcga = readtable('..\Data\tcga_median_tpm.txt');
arrayData.genes = tcga.gene;
arrayData.tissues = cell({ 'CRC'}); %regexprep(tcga.Properties.VariableNames(2:end)','_',' ');
arrayData.levels = median(table2array(tcga(:,{ 'COADTP', 'READTP'})),2);

load('..\Model\Human-GEM.mat')

model = ihuman;
model = addBoundaryMets(model);

taskStruct = parseTaskList('..\Data\metabolicTasks_Essential.txt');
%blockedReactions = findBlockedReaction(model);
[~,deletedDeadEndRxns] = simplifyModel(model,true,false,true,true,true);

%cModel = removeReactions(model,deletedDeadEndRxns,false,true);
%c

% add pre-processing results to arrayData structure
arrayData.deletedDeadEndRxns = deletedDeadEndRxns;
%arrayData.taskReport = taskReport;
%arrayData.essentialRxnMat = essentialRxnMat;
arrayData.threshold = 1;

%% Run tINIT algorithm

% initialize params input
params = {};

n_models = numel(arrayData.tissues);

% initialize INIT output structure
tissue_model = {};
j = 1;

% First try to run tINIT with shorter time limit. If it fails, then
% try again with a longer time limit.
disp(arrayData.tissues{j})
try
    params.TimeLimit = 1000;
    init_model = getINITModel2(model,arrayData.tissues{j},[],[],arrayData,[],true,[],true,true,taskStruct,params);
catch
    params.TimeLimit = 5000;
    init_model = getINITModel2(model,arrayData.tissues{j},[],[],arrayData,[],true,[],true,true,taskStruct,params);
end
init_model.b = zeros(length(init_model.b),1);
%init_model = removeDeadEnds(init_model);
% add ATP phosphohydrolase
%init_model = addReaction(init_model,'MAR03964','reactionName' ,'ATP phosphohydrolase','reactionFormula','ATP + H2O => ADP + H+ + Pi', 'subSystem','Transport reactions','checkDuplicate', 1);
if isempty(init_model.rxns(contains(init_model.rxns, 'MAR03964')))
    init_model = addReaction(init_model,'MAR03964','reactionName' ,'ATP phosphohydrolase','reactionFormula','MAM01371c + MAM02040c => MAM01285c + MAM02039c + MAM02751c', 'subSystem','Transport reactions','checkDuplicate', 1);
end
tissue_model.tissue_full_model{j,1} = init_model;
%simplify reaction
init_model = simplifyModel(init_model);
init_model = identifyBlockedRxns(init_model);

init_model.id = arrayData.tissues{j};

% save results
save('..\Model\tinit_model','init_model');
