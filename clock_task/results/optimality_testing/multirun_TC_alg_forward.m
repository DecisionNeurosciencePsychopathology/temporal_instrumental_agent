function [totcost, costs, seeds] = multirun_TC_alg_forward(conds, runseed, args)
    %use multiple runs of data to identify optimal parameters for TC model.
    rng(runseed);
    seeds=randi([1 500], length(conds), 1); %just the seed for reward outcomes
    
    costs=zeros(length(conds), 1);
    %execute runs in parallel
    parfor i = 1:length(conds)
        thiscall=args; %need to make a local version of args for parpool to work
        thiscall{4} = seeds(i,:); %rng seed is 4th argument
        thiscall{3} = conds{i}; %third argument is contingency to run
        costs(i) = TC_Alg_forward(thiscall{:});
    end
    
    totcost=sum(costs);
end

%function [cost, RTpred, ret]=TC_Alg_forward(params, priors, cond, rngseeds, ntrials, rtbounds)