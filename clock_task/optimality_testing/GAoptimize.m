function [Xga,Fga] = GAoptimize(agent, runarray)

%optionsGA = gaoptimset('PlotFcns',@gaplotbestf,'Display','iter','InitialPopulation',agent.init_params,...
%    'UseParallel','always','Vectorized', 'off','TolFun',1e-5);

optionsGA = gaoptimset('InitialPopulation',agent.init_params, ...
    'UseParallel',false,'Vectorized', 'off','TolFun',1e-4); %, ...
%    'PlotFcns',@gaplotbestf,'Display','iter');

%to pass additional parameters to fitness function using ga, an anonymous function is suggested
%other parameters are passed in based on current value at the time of function definition
f = @(x)fitnessfcn_optimize(x, agent, runarray);

try
    debug=agent.debug;
catch
    debug=0; %debug off by default
end

if debug
    %put here to to a fast test of parfor (just return initial pars)
    Xga=agent.init_params;
    Fga=rand;
else
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
end
