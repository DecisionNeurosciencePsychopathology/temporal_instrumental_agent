function [Xga,Fga] = GAoptimize(init_params,lower_bounds,upper_bounds,numVars,agt,cnd)

% numVars = 2; %Just epsilon(c) and the spread
% init_params_1 = [0.6    0.08 ];
% lower_bounds = [0.1 .01];
% upper_bounds = [.99 0.25];
optionsGA = gaoptimset('PlotFcns',@gaplotbestf,'Display','iter','InitialPopulation',init_params,...
    'UseParallel','always','Vectorized', 'off','TolFun',1e-5);

global agentNOW
global conditionNOW
agentNOW = agt;
conditionNOW = cnd;
[Xga,Fga] = ga(@fitnessfcn_optimize, ...
numVars, ... %number of free parameters
[], ... A*x
[], ... %b
[], ... %Aeq
[], ... %beq
lower_bounds, ... %LB
upper_bounds, ... %UB
[], ... %nonlcon
optionsGA);
