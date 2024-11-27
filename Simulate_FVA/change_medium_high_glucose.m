function [model,basisMedium]  = change_medium_high_glucose(model, closed)
%changeMedia_batch
% https://www.thermofisher.com/uk/en/home/technical-resources/media-formulation.9.html

if nargin < 2
    closed = true;
end

% 11995 - DMEM, high glucose, pyruvate
medium_composition={
'MAR09067'
'MAR09066'
'MAR09065'
'MAR09063'
'MAR09038'
'MAR09039'
'MAR09040'
'MAR09041'
'MAR09042'
'MAR09043'
'MAR09069'
'MAR09044'
'MAR09045'
'MAR09064'
'MAR09046'
'MAR09083'
'MAR09145'
'MAR09146'
'MAR09378'
'MAR09143'
'MAR09159'
'MAR09361'
'MAR09082'
'MAR09096'
'MAR09074'
'MAR09077'
'MAR09072'
'MAR09034'
'MAR09133'


};

% Medium concentrations
met_Conc_mM=[
0.4
0.39810428
0.20127796
4
0.2
0.8015267
0.8015267
0.7978142
0.20134228
0.4
0.4
0.79831934
0.078431375
0.39846742
0.8034188
0.028571429
0.008385744
0.009070295
0.032786883
0.00106383
0.011869436
0.04
1.8018018
0.000247525
0.8139166
110.344826
0.9057971
25
1


];


current_inf = 1000;
set_inf =1000;

cellConc = 2.17*1e6;
t= 24;
cellWeight = 3.645e-12;


%% Definition of basic medium (defines uptake from the medium, not captured by the medium composition, all with same constraints)
mediumCompounds = {'MAR09058'; 'MAR09079'; 'MAR09047'; 'MAR09078'; 'MAR11420'; 'MAR09048'; 'MAR09072'; 'MAR09074' };

mediumCompounds_lb = -1000;
%mediumCompounds = {};
%mediumCompounds_lb = -1000;

%% Determine constraints to simulate medium: 
[model,basisMedium] = setMediumConstraints(model, set_inf, current_inf, medium_composition, met_Conc_mM, cellConc, t, cellWeight, mediumCompounds, mediumCompounds_lb);


if closed
    %constrain consumption of exchange metabolites to 0.
    compartment = 'e';
    compartmentMets = ~cellfun(@isempty, strfind(model.mets, compartment));
    % Find reactions that involve the above mets
    compartmentRxns = model.rxns(any(model.S(compartmentMets, :)));
    compartmentReactions = [compartmentRxns, printRxnFormula(model, 'rxnAbbrList', compartmentRxns, 'printFlag', false)];
    
    compartmentReactions(ismember(compartmentReactions,medium_composition))='';
    compartmentReactions(ismember(compartmentReactions,mediumCompounds))='';

    model.lb(find(ismember(model.rxns,compartmentReactions)))= 0; 
    %for i = 1:length(compartmentReactions)
    %    model.lb(contains(model.rxns, compartmentReactions{i,1})) = 0;
    %    %disp(model.lb(contains(model.rxns, compartmentReactions{i,1})))
    %end
end

end