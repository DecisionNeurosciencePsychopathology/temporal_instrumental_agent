function [Xga,Fga] = GAoptimize(agent, runarray)

% optionsGA = gaoptimset('PlotFcns',@gaplotbestf,'Display','iter','InitialPopulation',init_params,...
%     'UseParallel','always','Vectorized', 'off','TolFun',1e-5);

optionsGA = gaoptimset('InitialPopulation',agent.init_params, ...
		       'UseParallel',false,'Vectorized', 'off','TolFun',1e-5); %, ...
    %'PlotFcns',@gaplotbestf,'Display','iter');

%to pass additional parameters to fitness function using ga, an anonymous function is suggested
%other parameters are passed in based on current value at the time of function definition
f = @(x)fitnessfcn_optimize(x, agent, runarray);

[Xga,Fga] = ga(f, ...
agent.k, ... %number of free parameters
[], ... A*x
[], ... %b
[], ... %Aeq
[], ... %beq
agent.lower_bounds, ... %LB
agent.upper_bounds, ... %UB
[], ... %nonlcon
optionsGA);

%put here to to a fast test of parfor
%Xga=agent.init_params;
%Fga=rand;
