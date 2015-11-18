function [totcost, costs, seeds] = multirun_cliff_walker_optimize(params, agent, runarray)
%use multiple runs of data to identify optimal parameters for Q-learning or SARSA agent
nruns=length(runarray);

costs=NaN(nruns, 1);
%execute runs in parallel
for i = 1:nruns
    thiscall={agent.clock_options, runarray(i), runarray(i).seeds, params}; %need to make a local version of args for parpool to work
    costs(i) = ClockWalking_3D_discountedEv_optimize(thiscall{:});
end

totcost=sum(costs);
end