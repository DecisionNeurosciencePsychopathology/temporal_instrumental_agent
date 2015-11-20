function [totcost, costs, seeds] = multirun_clock_logistic_operator(nruns, runseed, args)
    %use multiple runs of data to identify optimal parameters for logistic
    %operator.
    rng(runseed);
    seeds=randi([1 500], nruns, 2);
    
    costs=zeros(nruns, 1);
    %execute runs in parallel
    parfor i = 1:nruns
        thiscall=args; %need to make a local version of args for parpool to work
        thiscall{2} = seeds(i,:); %second argument is seed
        costs(i) = clock_logistic_operator(thiscall{:});
    end
    
    totcost=sum(costs);
end