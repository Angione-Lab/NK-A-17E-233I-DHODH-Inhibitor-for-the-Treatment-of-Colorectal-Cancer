clear all
clc

[num,txt,~] = xlsread("../Data/genes_fpkm.csv");
genes_in_dataset = txt(2:end,1);
conditions = txt(1,2:end);
gene_exp = num(:, 1:end);

load('..\Model\Human-GEM.mat')
model = ihuman;
[model, removedMets, removedRxns] = removeDeadEnds(model);
model.c(contains(model.rxns,'MAR13082')) = 1;


%Dehydroorotate  =>MAM03103c
%Orotate         => MAM02659c
%Glutamine       => MAM01975c
%Pyruvate        => MAM02819c
%Citrate         =>MAM01587c

%NAD             => MAM02552c
%NADH            =>MAM02553c
%Aspartate       =>MAM01370c
%Acetyl CoA      =>MAM01261c

%ubiquinol          =>MAM03102c
%ubiquinone         =>MAM03103c

% Add sink reactions
model = addSinkReactions(model, {'MAM03103c', 'MAM02659c', 'MAM01975c','MAM02819c', 'MAM01587c', 'MAM02552c', 'MAM02553c', 'MAM01370c', 'MAM01261c', 'MAM03102c', 'MAM03103c'});

% set medium constrain
[model,basisMedium]  = change_medium_high_glucose(model);

genes_in_model = model.genes;
gene_exp = gene_exp(contains(genes_in_dataset, genes_in_model),:);
genes_in_dataset = genes_in_dataset(contains(genes_in_dataset, genes_in_model));

% define parameters over which to iterate
gamma = [2];
threshold = [25];

[reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model);

genes = model.genes;
Cmin{size(gene_exp, 2), length(model.rxns)} = [];
Cmax{size(gene_exp, 2), length(model.rxns)} = [];
GeneExpressionArray = ones(numel(genes),1); 

k = 1;
tic

for g = 1:numel(gamma)                
    changeCobraSolver('gurobi', 'all');
    %changeCobraSolverParams('QP', 'feasTol', 1e-3);
    %changeCobraSolverParams('LP', 'feasTol', 1e-3);
    %changeCobraSolverParams('QP', 'method', 1);
    %changeCobraSolver('ibm_cplex', 'all');
    
	gam = gamma(g);
    new_k = main(threshold, gam, gene_exp,genes_in_dataset,conditions,model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k);
    k = new_k;
end
toc

save(['Results\CRC'])

function new_k = main(threshold, gam, gene_exp,genes_in_dataset,conditions,t_model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k)
    gamma = gam;
    if isempty(gcp('nocreate'))
        parpool();   
    end
    
    pool = gcp('nocreate'); % check if pool was successfully created
    if isempty(pool) % if there is no active pool, then throw an error
        error('\nError! No parallel pool created!')
    end
    for tr = 1:numel(threshold)                                                        
        data = gene_exp ;
       % applying the bounds
        fprintf('Iteration (k): %d, Gamma: %d, Threshold: %d\n',k,gamma(g),threshold(tr));
        for t=1:size(data,2)          % in here we select a unique profile(one patient at the time)
           	disp(t)
            model = t_model;
            % set biomass as objective function
            model = setParam(model, 'obj', 'MAR13082', 1);

            expr_profile = data(:,t);
            expr_profileCorrected = prctile(expr_profile,99);
            expr_profile(expr_profile > expr_profileCorrected) = expr_profileCorrected;
            expr_profile = expr_profile./max(expr_profile); 
            
            pos_genes_in_dataset = zeros(numel(genes),1);% gene in the model human 1
            for i=1:length(pos_genes_in_dataset)
                position = find(strcmp(genes{i},genes_in_dataset),1); 
                if ~isempty(position)                                   
                    pos_genes_in_dataset(i) = position(1);              
                    GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));         
                end
            end
            if or(sum(isinf(GeneExpressionArray)) >= 1, sum(isnan(GeneExpressionArray)) >= 1)
                fprintf('\nError in the gene expression data!');
            end
            [minfluxes, maxfluxes] = transcriptomic_bounds(conditions(t), gamma(g), GeneExpressionArray, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length);
            Cmin(k,:) = [conditions(k), num2cell(transpose(minfluxes))];
            Cmax(k,:) = [conditions(k), num2cell(transpose(maxfluxes))];
            k = k +1;
        end
    end
    new_k = k; 
    Tmin = cell2table(Cmin);
    Tmin.Properties.VariableNames = vertcat({'conditions'}, model.rxns);

    Tmax = cell2table(Cmax);
    Tmax.Properties.VariableNames = vertcat({'conditions'}, model.rxns); 
    writetable(Tmin,'Results\FVA_min.csv');     
    writetable(Tmax,'Results\FVA_max.csv');     % set the folder where you want the data to be saved
end