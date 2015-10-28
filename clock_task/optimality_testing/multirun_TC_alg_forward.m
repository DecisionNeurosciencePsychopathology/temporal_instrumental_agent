function [totcost, costs, seeds] = multirun_TC_alg_forward(params, agent, runarray)
    %use multiple runs of data to identify optimal parameters for TC model.
    rng(agent.runseed);
    nruns=length(runarray);
    seeds=randi([1 500], nruns, 1); %just the seed for reward outcomes
    
    rtbounds=[0 5000]; %RT space for testing
    costs=NaN(nruns, 1);
    priors.V = 0; %don't give agent any insight into previous values (no SCEPTIC receives this boost)
    priors.Go = 0;
    prior.NoGo = 0;
    %execute runs in parallel
    for i = 1:nruns
        thiscall = {params, priors, runarray(i), seeds(i,:), agent.ntrials, rtbounds};
        costs(i) = TC_Alg_forward(thiscall{:});
    end
    
    totcost=sum(costs);
end