

clear all
clc

fluxFoldChange_pathways_all = [];
load ('Results\CRC.mat');

tbl_maxFlux = rows2vars(readtable('Results/FVA_max.csv'));
maxFlux = cell2mat(table2cell(tbl_maxFlux(2:end, 2:end)));

tbl_minFlux = rows2vars(readtable('Results/FVA_min.csv'));
minFlux = cell2mat(table2cell(tbl_minFlux(2:end, 2:end)));

%%

%%
%reactions in upregulated and downregulated genes

%% get reactions containing upregulated and downregulated genes
[model] = generateRules(model);
[results, ListResults] = findRxnsFromGenes(model, {'ENSG00000102967'}, 1,1);
%ListResults = ListResults(:, 1:5);
%
%%
ferroptosis_gene_list = {'ENSG00000138449', 'ENSG00000168003', 'ENSG00000001084', 'ENSG00000023909', 'ENSG00000197142', 'ENSG00000100983', 'ENSG00000151012'};
[results, ListResults] = findRxnsFromGenes(model, ferroptosis_gene_list, 1,1);
%% get significant extracellular transport reactions containing upregulated and downregulated genes
%ListResults = ListResults(ismember(ListResults(:,1), metRxns), :);
%ListResults = ListResults(strcmp(string(ListResults(:,3)), 'Transport reactions'), :);


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
rxnGene = [{'Reaction_id', 'Reaction Names', 'Pathways', 'Genes'}; rxn_genes]; 
xlswrite(['Results\Reactions_genes'], rxnGene)

%% Calculate reaction activity fold changes
% set to 0 values below cplex feasibility tolerance
maxFlux(maxFlux < 1e-06) = 0;
minFlux(minFlux < 1e-06) = 0;
% split reversible reactions
revIDs = model.lb < 0;
rxns = [model.rxns; strcat(model.rxns(revIDs), '_b')];
rxnNames = [model.rxnNames; strcat(model.rxnNames(revIDs), ' (backward)')];
subSystems = [model.subSystems; model.subSystems(revIDs)];
formulas = printRxnFormula(model, model.rxns, false);
formulas = strrep(formulas, '[]', '');
formulas = [formulas;formulas(revIDs)];

gprRule = model.grRules;
gprRule = strrep(gprRule, '[]', '');
gprRule = [gprRule; gprRule(revIDs)];

maxFlux = [maxFlux; -minFlux(revIDs, :)];
maxFlux(maxFlux < 0) = 0;
maxFlux(sum(diff(sort(maxFlux,2),1,2)~=0,2)+1 <=3,1:end) = 0;


% re-sort reactions
[rxns, sort_idx] = sort(rxns);
rxnNames = rxnNames(sort_idx);
subSystems = subSystems(sort_idx);
formulas = formulas(sort_idx);
gprRule = gprRule(sort_idx);
maxFlux = maxFlux(sort_idx, :);
correctedMaxFlux = maxFlux + 1;

subSystems(1129) = {'Heme synthesis'};
subSystems = cellstr(string(subSystems));

% % get solution span for each flux
fluxFoldChange = mean(correctedMaxFlux(:, 4:6),2)./ mean(correctedMaxFlux(:, 1:3),2)  ;

fluxFoldChange(isnan(fluxFoldChange)) = 1;

metabolites = findMetsFromRxns(model, rxns);
compartments = getCompartment(metabolites);

%% Export reaction data

t1 = [{'Reaction ID', 'Reaction name', 'Formula', 'Pathway', 'GPR rules', 'ctrl_1', 'ctrl_2', 'ctrl_3', 'TTT_1', 'TTT_2', 'TTT_3','Average Maxflux wildtype','Average Maxflux m4', 'Foldchange m4_vs_wildtype', 'log2(Foldchange m4_vs_wildtype)'}; ...
    rxns, rxnNames, formulas, subSystems, gprRule, num2cell(correctedMaxFlux) ,num2cell(mean(correctedMaxFlux(:, 1:3),2)), num2cell(mean(correctedMaxFlux(:, 4:6),2)), num2cell(fluxFoldChange), num2cell(log2(fluxFoldChange))];


xlswrite(['Results/CRC_results'], t1, 'reactions');
%

%% Reactions enrichment
%rxns = cellfun(@(x,y) [x  ' (' y ') '],T{:,'rxns'},T{:,'comp'},'un',0) ;
%rxnNames = cellfun(@(x,y) [x  ' (' y ') '],T{:,'rxnNames'},T{:,'comp'},'un',0) ;
logfluxFoldChange_rxns = log2(fluxFoldChange);
upregulated = and(logfluxFoldChange_rxns > prctile(logfluxFoldChange_rxns, 95), logfluxFoldChange_rxns > 0.01);
downregulated = and(logfluxFoldChange_rxns <= prctile(logfluxFoldChange_rxns, 5), logfluxFoldChange_rxns < -0.01);

up_rxns = rxns(upregulated(:,1));
down_rxns = rxns(downregulated(:,1));
uniquerxns = unique([up_rxns; down_rxns]); %reactions
enriched_subsystems = subSystems(ismember(rxns, uniquerxns));
uniquerxnNames = rxnNames(ismember(rxns,uniquerxns));
uniquerxnNames = strcat(uniquerxnNames, '(', uniquerxns , ')');
unique_gprRule = gprRule(ismember(rxns,uniquerxns));


fluxFoldChange_rxns = logfluxFoldChange_rxns(ismember(rxns,uniquerxns),1) ;

%sort by first FC value
[~, sort_idx] = sort(fluxFoldChange_rxns(:,1));
uniquerxns = uniquerxns(sort_idx, :);
uniquerxnNames = uniquerxnNames(sort_idx, :);
fluxFoldChange_rxns = fluxFoldChange_rxns(sort_idx, :);
enriched_subsystems = enriched_subsystems(sort_idx, :);
unique_gprRule = unique_gprRule(sort_idx, :);

%save reaction enrichment to file
t1 = [{'Reactions', 'Reaction Names','pathways', 'GPR', 'FC M4 vs wildtype'}; ...
    cell(uniquerxns),cell(uniquerxnNames) cell(enriched_subsystems),cell(unique_gprRule), num2cell(fluxFoldChange_rxns)];
xlswrite(['Results/reaction_foldchange'], t1, 'reactions');


%% Plot pathway fold changes
L = 10; % number of datapoints
indexValue = 0; % value for which to set a particular color
topColor = [0.5 0 0]; % color for maximum data value (dark red = [0.5 0 0])
indexColor1 = [1 0 0]; % color for intermediate data value (red = [1 0 0])
indexColor = [1 1 1]; % color for null value (white = [1 1 1])
bottomcolor = [0 0 1]; % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and maximum values
largest = max(max(fluxFoldChange_rxns));
smallest = min(min(fluxFoldChange_rxns));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multiplying number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),10*index)',...
            linspace(bottomcolor(2),indexColor(2),10*index)',...
            linspace(bottomcolor(3),indexColor(3),10*index)'];
customCMap3 = [linspace(indexColor(1),indexColor1(1),10*(index))',...
            linspace(indexColor(2),indexColor1(2),10*(index))',...
            linspace(indexColor(3),indexColor1(3),10*(index))'];
% Create color map ranging from index color to top color
% Multiplying number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor1(1),topColor(1),10*(L-2*index))',...
            linspace(indexColor1(2),topColor(2),10*(L-2*index))',...
            linspace(indexColor1(3),topColor(3),10*(L-2*index))'];
customCMap = [customCMap1; customCMap3; customCMap2];  % Combine colormaps

subplot('Position', [0.45 0.01 0.15 0.99])
h1 = heatmap({'FC m4 vs wildtype'}, uniquerxnNames, fluxFoldChange_rxns, 'Colormap', customCMap, ...
    'ColorLimits', [smallest largest]);
h1.NodeChildren(3).XAxis.TickLabelInterpreter = 'None';
h1.NodeChildren(3).YAxis.TickLabelInterpreter = 'None';
h1.NodeChildren(3).Title.Interpreter = 'None';

set(gca, 'FontSize', 3, 'FontName', 'Arial')
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperOrientation = 'portrait';
fig.PaperSize = [10, 17];
print('-fillpage', 'Results/crc_reactions_high_glucose', '-dpdf');


%% Pathway statistics
[fluxFoldChange_pathways, fluxFoldChange_pathways_sem, fluxFoldChange_pathways_std] = grpstats(fluxFoldChange, subSystems, {'mean', 'sem', 'std'});
%log2 transform mean fold changes
logfluxFoldChange_pathways = log2(fluxFoldChange_pathways);
uniqueSubSystems = unique(subSystems, 'stable');
[uniqueSubSystems, sort_idx] = sort(uniqueSubSystems);
fluxFoldChange_pathways = fluxFoldChange_pathways(sort_idx, :);
fluxFoldChange_pathways_sem = fluxFoldChange_pathways_sem(sort_idx, :);
fluxFoldChange_pathways_std = fluxFoldChange_pathways_std(sort_idx, :);
logfluxFoldChange_pathways = logfluxFoldChange_pathways(sort_idx, :);
t2 = [{'Pathway', 'Foldchange mean m4_vs_wildtype','Foldchange sem m4_vs_wildtype','Foldchange std m4_vs_wildtype', 'log2FC mean m4_vs_wildtype'}; ...
    uniqueSubSystems, num2cell(fluxFoldChange_pathways), num2cell(fluxFoldChange_pathways_sem), num2cell(fluxFoldChange_pathways_std),num2cell(logfluxFoldChange_pathways)];
xlswrite(['Results/CRC_results'], t2, 'pathways');

%% Reaction enrichment
a = and(log2(fluxFoldChange) > prctile(log2(fluxFoldChange), 95), log2(fluxFoldChange) > 0);
enrichment_table = FEA(a(:, 1), rxns);

a = and(log2(fluxFoldChange) <= prctile(log2(fluxFoldChange), 5), log2(fluxFoldChange) < 0);
enrichment_table = [enrichment_table; FEA(a(:, 1), rxns)];

enrichment_table = [repmat({''}, size(enrichment_table, 1), 1), enrichment_table];
enrichment_table(strcmp(enrichment_table(:, 2), 'p-value'), 1) = {'UP (FC > 95th percentile and > 0)'; 'DOWN (FC < 5th percentile and < 0)'};
xlswrite(['Results_sink/CRC_results'], enrichment_table, 'rxn_enrichment_table');


%% Pathway enrichment
a = and(log2(fluxFoldChange) > prctile(log2(fluxFoldChange), 85), log2(fluxFoldChange) > 0);
enrichment_table = FEA(a(:, 1), subSystems);

%
a = and(log2(fluxFoldChange) <= prctile(log2(fluxFoldChange), 15), log2(fluxFoldChange) < 0);
enrichment_table = [enrichment_table; FEA(a(:, 1), subSystems)];

enrichment_table = [repmat({''}, size(enrichment_table, 1), 1), enrichment_table];
enrichment_table(strcmp(enrichment_table(:, 2), 'p-value'), 1) = {'UP (FC > 85th percentile and > 0)'; 'DOWN (FC < 15th percentile and < 0)'};
xlswrite(['Results/CRC_results'], enrichment_table, 'enrichment_table');

% metabolite info to interpret the formulas
model.mets = strrep(model.mets, '[]', '');
t4 = [{'Metabolite ID', 'Metabolite name'}; ...
    model.mets, model.metNames];
xlswrite(['Results/CRC_results'], t4, 'metabolites');


% aggregate pathway-level mean fold changes
fluxFoldChange_pathways_all = [fluxFoldChange_pathways_all logfluxFoldChange_pathways];

%%
upregulated = and(logfluxFoldChange_pathways > prctile(logfluxFoldChange_pathways, 85), logfluxFoldChange_pathways > 0);
downregulated = and(logfluxFoldChange_pathways <= prctile(logfluxFoldChange_pathways, 15), logfluxFoldChange_pathways < 0);

up_pathways = uniqueSubSystems(upregulated(:,1));
down_pathways = uniqueSubSystems(downregulated(:,1));
uniqueSubSystems1 = unique([up_pathways; down_pathways]); %uniqueSubSystems
fluxFoldChange_pathways1 = logfluxFoldChange_pathways(ismember(uniqueSubSystems,uniqueSubSystems1),1);


%sort by first FC value
[~, sort_idx] = sort(fluxFoldChange_pathways1(:,1));
uniqueSubSystems1 = uniqueSubSystems1(sort_idx, :);
fluxFoldChange_pathways1 = fluxFoldChange_pathways1(sort_idx, :);


%% Plot pathway fold changes
L = 10; % number of datapoints
indexValue = 0; % value for which to set a particular color
topColor = [0.5 0 0]; % color for maximum data value (dark red = [0.5 0 0])
indexColor1 = [1 0 0]; % color for intermediate data value (red = [1 0 0])
indexColor = [1 1 1]; % color for null value (white = [1 1 1])
bottomcolor = [0 0 1]; % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and maximum values
largest = max(max(fluxFoldChange_pathways_all));
smallest = min(min(fluxFoldChange_pathways_all));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multiplying number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];
customCMap3 = [linspace(indexColor(1),indexColor1(1),100*(index))',...
            linspace(indexColor(2),indexColor1(2),100*(index))',...
            linspace(indexColor(3),indexColor1(3),100*(index))'];
% Create color map ranging from index color to top color
% Multiplying number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor1(1),topColor(1),100*(L-2*index))',...
            linspace(indexColor1(2),topColor(2),100*(L-2*index))',...
            linspace(indexColor1(3),topColor(3),100*(L-2*index))'];
customCMap = [customCMap1; customCMap3; customCMap2];  % Combine colormaps

% split pathways for the two heatmaps
idx = contains(uniqueSubSystems, 'translation and secretion', 'IgnoreCase', true);
uniqueSubSystems2 = uniqueSubSystems(idx);
fluxFoldChange_pathways2 = fluxFoldChange_pathways_all(idx, :);
% fix entries
% remove Protein translation and secretion
idx = contains(uniqueSubSystems2, 'Protein translation and secretion', 'IgnoreCase', true);
uniqueSubSystems2(idx) = [];
fluxFoldChange_pathways2(idx) = [];
for i = 1:length(uniqueSubSystems2)
    uniqueSubSystems2{i} = strrep(uniqueSubSystems2{i}, ' translation and secretion', '');
end
uniqueSubSystems2 =  uniqueSubSystems2(fluxFoldChange_pathways2 ~= 0);
fluxFoldChange_pathways2 = fluxFoldChange_pathways2(fluxFoldChange_pathways2 ~= 0);
% re-sort entries
[uniqueSubSystems2, sort_idx] = sort(uniqueSubSystems2);
fluxFoldChange_pathways2 = fluxFoldChange_pathways2(sort_idx, :);

subplot('Position', [0.35 0.01 0.15 0.99])
h1 = heatmap({'FC m4 vs wildtype'}, uniqueSubSystems1, fluxFoldChange_pathways1, 'Colormap', customCMap, ...
    'ColorLimits', [smallest largest]);
set(gca, 'FontSize', 3, 'FontName', 'Arial')
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperOrientation = 'portrait';
fig.PaperSize = [10, 17];
print('-fillpage', 'Results/CRC_results', '-dpdf');

%%
function resultCell = FEA(rxnSet, group)

if nargin < 2
    error('The function FEA must be called with reaction set and group as arguments')
end
if ~isvector(rxnSet)
    error('Please provide the indices of the reactions e.g. 1:10')
end
if ~iscell(group)
    error('Please provide the group name as cell array of characters e.g. the subSystem field of any metabolic model ')
end

%Temporary Warning until FEA statistics are checked for multiple classes.
if iscell(group{1}) %Potentially multiple subSystems
    if any(cellfun(@numel, group) > 2)
        warning('Multiple subSystems detected for some reactions. FEA statistics might not be correct.\n Please consider using only one subSystem per reaction.')
    end
end

% compute frequency of enriched terms
%groups = eval(['model.' group]);
groups = group;
if iscell([groups{:}])
   [uniquehSubsystemsA] = unique([groups{:}]);
   presenceindicator = false(numel(uniquehSubsystemsA),numel(group));
   for i = 1:numel(groups)
       presenceindicator(:,i) = ismember(uniquehSubsystemsA,groups{i});
   end   
   [K,~] = find(presenceindicator);
else
    %This works only for fields which have a single entry.
    [uniquehSubsystemsA, ~, K] = unique(groups);
end
% fetch group
%enRxns = eval(['model.' group '(rxnSet)']);
enRxns = group(rxnSet);
m = length(uniquehSubsystemsA);
allSubsystems = zeros(1, m);

% look for unique occurences
if iscell([enRxns{:}])
   [uniquehSubsystems] = unique([enRxns{:}]);
   presenceindicator = false(numel(uniquehSubsystems),numel(group));
   for i = 1:numel(enRxns)       
        presenceindicator(:,i) = ismember(uniquehSubsystems,enRxns{i});
   end   
   [J,~] = find(presenceindicator);
else
    %This works only for fields which have a single entry.
    [uniquehSubsystems, ~, J] = unique(enRxns);
end

occ = histc(J, 1:numel(uniquehSubsystems));
[l, p] = intersect(uniquehSubsystemsA, uniquehSubsystems);
allSubsystems(p) = occ;

% compute total number of reactions per group
nRxns = histc(K, 1:numel(uniquehSubsystemsA));  % the number of reactions per susbsystem

% Compute p-values
% gopvalues = hygepdf(allSubsystems', max(nRxns), max(allSubsystems), nRxns);
gopvalues = hygecdf(allSubsystems'-1, repmat(sum(nRxns), length(nRxns), 1), nRxns, repmat(sum(allSubsystems), length(nRxns), 1), 'upper');

% take out the zeros for one-sided test
nonZerInd = find(allSubsystems);

% sort p-values
[m, rxnInd] = sort(gopvalues);

% intersect non zero sets with ordered pvalues
[~, nonZeroInd] = intersect(rxnInd, nonZerInd);
orderedPval = rxnInd(sort(nonZeroInd));

% Build result cell
% initilize variable
resultCell = cell(length(orderedPval) + 1, 5);
resultCell(1, :) = {'p-value', 'Adjusted p-value', 'Pathway', 'Enriched set size', 'Total set size'};

% P values
resultCell(2:end, 1) = num2cell(gopvalues(orderedPval));

% correct for multiple testing with FDR
resultCell(2:end, 2) = num2cell(mafdr(cell2mat(resultCell(2:end, 1)), 'BHFDR', true));

% Group name
resultCell(2:end, 3) = uniquehSubsystemsA(orderedPval);

% Test size
resultCell(2:end, 4) = num2cell(allSubsystems(orderedPval))';

% Total group size
resultCell(2:end, 5) = num2cell(nRxns(orderedPval));

end

