function [totcost, costs, seeds] = multirun_TC_alg_forward(params, agent, runarray)
    %use multiple runs of data to identify optimal parameters for TC model.
    nruns=length(runarray);
    
    rtbounds=[0 5000]; %RT space for testing
    costs=NaN(nruns, 1);
    priors.V = 0; %don't give agent any insight into previous values (no SCEPTIC receives this boost)
    priors.Go = 0;
    priors.NoGo = 0;
    %execute runs in parallel
    for i = 1:nruns
        thiscall = {params, priors, runarray(i), runarray(i).seeds, agent.ntrials, rtbounds};
        costs(i) = TC_Alg_forward(thiscall{:});
    end
    
    totcost=sum(costs);
end
