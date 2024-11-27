

clear all
clc

fluxFoldChange_pathways_all = [];
load ('Results_fpkm\CRC.mat');

tbl_maxFlux = rows2vars(readtable('Results/FVA_max.csv'));
maxFlux = cell2mat(table2cell(tbl_maxFlux(2:end, 2:end)));

tbl_minFlux = rows2vars(readtable('Results/FVA_min.csv'));
minFlux = cell2mat(table2cell(tbl_minFlux(2:end, 2:end)));

%% get reactions containing upregulated and downregulated genes
[model] = generateRules(model);
[results, ListResults] = findRxnsFromGenes(model, {'ENSG00000102967'}, 1,1);
%ListResults = ListResults(:, 1:5);

%% update reaction names
rxn_file = readtable('../Data/reactions.tsv.txt');
rxn_file = table2cell(rxn_file);
model.reactions = model.rxnNames;
for c = 1:height(model.rxns)
    if cellfun(@isempty,model.rxnNames(c))
        %model.rxnNames(c) = model.rxns(c); 
        reaction = model.rxns(c);
        reconid = rxn_file(contains(rxn_file(:,1), reaction),7);
        if ~cellfun(@isempty,reconid) 
            disp(reconid)
            model.rxnNames(contains(model.rxns, reaction)) = reconid;
        end
    end 
end
for c = 1:height(model.rxns)
    if cellfun(@isempty,model.rxnNames(c))
       model.rxnNames(c) = model.rxns(c); 
    end 
end

%% exchange metabolites
isExRxn = findExcRxns(model);
exchangeRxns = model.rxns(isExRxn);
max_exflux = maxFlux(isExRxn,:);
min_exflux = minFlux(isExRxn,:);
netExFlux = max_exflux +  min_exflux;

rxn_genes{length(model.rxns), 4} = [];
for i = 1:length(model.rxns)
    genes = findGenesFromRxns(model, model.rxns{i});
    rxn_genes{i,1} = model.rxns{i};
    rxn_genes{i,2} = model.rxnNames{i};
    rxn_genes{i,3} = model.subSystems{i}{1};
    rxn_genes{i,4} = strjoin(genes{1}, ', ');
end
%rxnGene = [{'Reaction_id', 'Reaction Names', 'Pathways', 'Genes'}; rxn_genes]; 
%xlswrite(['Results_fpkm\Reactions_genes'], rxnGene)

%% Calculate reaction activity fold changes
% set to 0 values below cplex feasibility tolerance
maxFlux(maxFlux < 1e-06) = 0;
minFlux(minFlux < 1e-06) = 0;
% split reversible reactions
revIDs = model.lb < 0;
rxns = model.rxns; %[model.rxns; strcat(model.rxns(revIDs), '_b')];
rxnNames = model.rxnNames; %[model.rxnNames; strcat(model.rxnNames(revIDs), ' (backward)')];
subSystems = model.subSystems;  %[model.subSystems; model.subSystems(revIDs)];
formulas = printRxnFormula(model, model.rxns, false);

gprRule = model.grRules;

maxFlux(maxFlux < 0) = 0;

correctedMaxFlux = maxFlux + 1;

% % get solution span for each flux
fluxFoldChange = mean(correctedMaxFlux(:, 4:6),2)./ mean(correctedMaxFlux(:, 1:3),2)  ;

fluxFoldChange(isnan(fluxFoldChange)) = 1;

metabolites = findMetsFromRxns(model, rxns);
compartments = getCompartment(metabolites);

%MAM00180c + MAM03103m  <=> MAM02659c + MAM03102m 
TTT_maxflux = mean(correctedMaxFlux(:, 4:6),2);
Ctrl_maxflux = mean(correctedMaxFlux(:, 1:3),2);

metAbbr = "MAM03102m";
N = 20;
metaboliteMassBalance(model, metAbbr, TTT_maxflux, N)

metaboliteMassBalance(model, metAbbr, Ctrl_maxflux, N)

