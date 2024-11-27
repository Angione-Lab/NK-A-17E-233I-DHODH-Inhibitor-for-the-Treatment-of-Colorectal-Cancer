function [minFlux, maxfluxes] = transcriptomic_bounds(condition, gamma, gene_expression, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length)
% function modified from evaluate_objective.m
yt=gene_expression';      
eval_reaction_expression = reaction_expression;

for i=ixs_geni_sorted_by_length 
    posizioni_gene = pos_genes_in_react_expr{i};
    for j=1:length(posizioni_gene)  
        eval_reaction_expression{posizioni_gene(j)} = strrep(eval_reaction_expression{posizioni_gene(j)}, genes{i}, num2str(yt(i),'%.15f'));  
    end
end

eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1.0'};  
num_reaction_expression = zeros(1,length(eval_reaction_expression));

for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};
    num_parenthesis = numel(strfind(str,')'));
    while (num_parenthesis > 32) 
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
        substrings_to_replace = regexp(str, to_replace, 'match');
        if isempty(substrings_to_replace)
            num_parenthesis = 0; 
        else
            for j = 1:numel(substrings_to_replace)
                ss_rep = substrings_to_replace{j};
                str = strrep(str,ss_rep,num2str(eval(ss_rep),'%.15f'));
            end
            num_parenthesis = numel(strfind(str,')'));
        end
    end
    str = regexprep(str,'/','');
    try
        num_reaction_expression(i) = eval(str);
    catch
        num_reaction_expression(i) = 1;
    end
end

%reaction_expressionUpCorrected = prctile(num_reaction_expression,99);
%num_reaction_expression(num_reaction_expression > reaction_expressionUpCorrected) = reaction_expressionUpCorrected;

if or(sum(isinf(num_reaction_expression)) >= 1, sum(isnan(num_reaction_expression)) >= 1)
    fprintf('\nError in the evaluated and corrected gene expression data!');
end

for i=1:length(num_reaction_expression)   %loop over the array of the geneset expressions
    if ~isnan(num_reaction_expression(i))
        if num_reaction_expression(i) < 1
            model.lb(i) = model.lb(i)*(gamma * num_reaction_expression(i)); % model default x gene expression to the power of gamma
            model.ub(i) = model.ub(i)*(gamma * num_reaction_expression(i));
        end
    end
end

if or(sum(isinf(model.lb)) >= 1, sum(isnan(model.lb)) >= 1)
    fprintf('\nError in the new lower bounds!');
end
if or(sum(isinf(model.ub)) >= 1, sum(isnan(model.ub)) >= 1)
    fprintf('\nError in the new upper bounds!');
end

%% Compute the fluxes
start = tic;
[minFlux, maxfluxes] = fluxVariability(model, 100, 'max', model.rxns, 1);

finish = toc(start);
fprintf('\nTime taken: %f minutes', finish/60);

format longG; format compact;
