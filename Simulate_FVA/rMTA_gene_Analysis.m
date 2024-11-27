%clear
%initCobraToolbox(false) % false, as we don't want to update
transcript_separator = '.';

changeCobraSolver('ibm_cplex', 'all');
%changeCobraSolver('gurobi', 'all');
%load('..\Model\Human-GEM.mat')
%model = ihuman;

% we will load tinit model to generate flux samples
load('..\Model\tinit_model.mat')
model = init_model;
model.c(contains(model.rxns,'MAR13082')) = 1;
%[model, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

% set medium constrain
cd ..\Simulate_FVA
[model,basisMedium]  = change_medium_high_glucose(model);

cd ..\rMTA

%%
% Sampling Method and options
sampling_method = 'ACHR';	%{('CHRR'), 'ACHR'}
sampling_options = struct();
sampling_options.nWarmupPoints = 5000;      	% (default)
sampling_options.nPointsReturned  = 2000;   	% (default)
sampling_options.nStepsPerPoint  = 500;    	% (default = 200)

% Now COBRAtoolbox include sampler
tic
if isempty(gcp('nocreate'))
    parpool();   
end

pool = gcp('nocreate'); % check if pool was successfully created
if isempty(pool) % if there is no active pool, then throw an error
    error('\nError! No parallel pool created!')
end

[modelSampling,samples] = sampleCbModel(model,'sampleFiles',sampling_method,sampling_options);
%load(['Data' filesep 'Sampling_results.mat']);
TIME.sampling = toc;

sampleStats = calcSampleStats(samples);
% from reduced index to model.rxns index
idx = zeros(size(samples,1),1);
for i = 1:numel(idx)
    try
        idx(i) = find(strcmp(model.rxns,modelSampling.rxns{i}));
    catch
        idx(i) = find(cellfun(@length,strfind(model.rxns,modelSampling.rxns{i})));
        sampleStats.mean(i)=-1*sampleStats.mean(i); %Those reactions are reversed;
    end
end
rxnInactive = setdiff(1:length(model.rxns),idx); % inactive reactions
fields = fieldnames(sampleStats);
for i = 1:numel(fields)
    aux = sampleStats.(fields{i});
    sampleStats.(fields{i}) = zeros(size(model.rxns));
    sampleStats.(fields{i})(idx) = aux;
    clear aux
end

aux = samples;
samples = zeros(size(model.rxns,1),sampling_options.nPointsReturned);
samples (idx,:) = aux;
clear aux;   

Vref = sampleStats.mean;

%%
% Differentially expressed genes
% Neccesary variables: 'gene','logFC','pval'
% 'Gene_ID' must be the same nomenclature as the metabolic model
% Study must be uploaded as DISEASE VS HEALTHY/CONTROL
filename_differentially_expressed_genes = 'DEG.csv';
logFC_requiered = 0; % change if necesary
pval_requiered = 0.1; % change if necesary
 
differ_genes = readtable(fullfile('..\Data',filename_differentially_expressed_genes),...
 'ReadVariableNames',true);
 
differ_genes.pval = differ_genes.padj; % requiered variable by code
differ_genes.gene = differ_genes.Var1; % requiered variable by code
 

% Here we should obtain an array similar to rxnHML, in which we have the
% information of whatever an expresion should increase, decrease or nothing
% (+1)R_f    (-1)R_b     (0)unchanged
% This vector is called rxnFBS (Forward, Backward, Unchanged)
 
rxnFBS = diffexprs2rxnFBS(model, differ_genes, Vref, ...
 'SeparateTranscript', transcript_separator, 'logFC', logFC_requiered, 'pval', pval_requiered);

%%
%% STEP 3: run rMTA algorithm, implemented in a function available in COBRA toolbox
% Both MTA and rMTA use the same parameters: epsilon and alpha:

% Define alpha values to calculate rMTA
alpha_values = [0.66];  % (default range of values) % better to have more values
%  It has been included 0.66 as it is the original value used in the
%  original paper ('Yizhak et al, 2013')
num_alphas = length(alpha_values);

% Calculate epsilon, different for each reaction and with a minimum required change of 1e-3 (%default)
epsilon = calculateEPSILON(samples, rxnFBS);
%% 
% One we have defined the parameters required by rMTA/MTA, we can use the 
% COBRA function to calculate the transformation score (TS score):

changeCobraSolver('ibm_cplex', 'all')
delete(pool)
if isempty(gcp('nocreate'))
    fprintf('Creating new pool')
    parpool();   
end
pool = gcp('nocreate'); % check if pool was successfully created
if isempty(pool) % if there is no active pool, then throw an error
    error('\nError! No parallel pool created!')
end
% execute the code
tic
[TSscore, deletedGenes, Vres] = rMTA(model, rxnFBS, Vref, alpha_values, epsilon, ...
    'timelimit', 60, 'SeparateTranscript', transcript_separator, 'printLevel', 1);
TIME.rMTA = toc ;
%% STEP 4: save results in an Excel for study
% First of all, we are going to clean the Results folder and save the results 
% from this tutorial

save(['CRC_case_genes.mat'])


%%
% gene information
filename3 = 'C:\Users\suraj\OneDrive - Teesside University\1.PhDworks\Ewing_GSMM\New_data_for_GSMM\GSMM\data\GeneInfo_HomoSapiens_ENSEMBL_103.txt';
biomart_genes = readtable(filename3,'ReadVariableNames',true);
biomart_genes.NCBIGene_formerlyEntrezgene_ID = cellfun(@num2str, num2cell(biomart_genes.NCBIGene_formerlyEntrezgene_ID), 'UniformOutput', 0);
[~, idx] = ismember(deletedGenes, biomart_genes.NCBIGene_formerlyEntrezgene_ID);
idx = idx (idx>0);
gene_info = biomart_genes(idx,:);
gene_info.gene = gene_info.NCBIGene_formerlyEntrezgene_ID;
geneID=table(deletedGenes, 'VariableNames', {'gene'});
gene_info = outerjoin(geneID,gene_info,'MergeKeys',true);

save(['CRC_simulation_rMTA_genes_results'])

rMTAsaveInExcel(['CRC__genes_case.xlsx'], TSscore, deletedGenes, alpha_values, ...
    'differ_genes', differ_genes)

 TIME
