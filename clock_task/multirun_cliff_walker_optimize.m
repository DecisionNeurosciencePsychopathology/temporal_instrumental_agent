function [totcost, costs, seeds] = multirun_cliff_walker_optimize(params, agent, runarray)
    %use multiple runs of data to identify optimal parameters for logistic
    %operator.
    rng(agent.runseed);
    nruns=length(runarray);
    seeds=randi([1 500], nruns, 1); %just the seed for reward outcomes
    
    costs=NaN(nruns, 1);
    %execute runs in parallel
    for i = 1:nruns        
        thiscall={agent.clock_options, runarray(i), seeds(i,:), params}; %need to make a local version of args for parpool to work
        costs(i) = ClockWalking_3D_discountedEv_optimize(thiscall{:});
    end
    
    totcost=sum(costs);
end


% [maxReward, constr, quits,cumReward,mov, ret] = ClockWalking_3D_discountedEv_optimize(options,m,rngseeds,params,plot_index,gra_options)