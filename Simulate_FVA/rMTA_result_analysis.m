clear all
clc

%top_genes = readtable('CRC_case.xlsx', Sheet='TOP_GENES_ID')

%top_rxns = readtable('CRC__rxn_case.xlsx', Sheet='TOP_GENES_ID');
%top_rxns = table2cell(top_rxns(:, 1))
[~, ~, top_rxns] = xlsread('CRC__rxn_case.xlsx', 'TOP_GENES_ID');
top_rxns = top_rxns(2:end,:);

load ('..\Simulate_FVA\Results\CRC.mat');

[model] = generateRules(model);

%[results, ListResults] = findRxnsFromGenes(model, table2cell(top_genes), 1,1);
%ListResults = cell2table(ListResults(:, 1:5));

formulas = printRxnFormula(model, top_rxns, false);
gprRule = model.grRules(contains(model.rxns, top_rxns));

a = [{'Pathways', 'Reactions', 'Reaction Names', 'Formula', 'grRules'}; ...
    model.subSystems(contains(model.rxns, top_rxns)), model.rxns(contains(model.rxns, top_rxns)),  model.rxnNames(contains(model.rxns, top_rxns)), formulas, gprRule, top_rxns]

writetable(cell2table(a), 'rmta_rxns_pathways.csv')


