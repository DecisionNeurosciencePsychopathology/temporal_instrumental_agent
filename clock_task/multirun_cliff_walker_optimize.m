function [totcost, costs, seeds] = multirun_cliff_walker_optimize(nruns, runseed, args, vperm)
    %use multiple runs of data to identify optimal parameters for logistic
    %operator.
    rng(runseed);
    seeds=randi([1 500], nruns, 2);
    
    costs=zeros(nruns, 1);
    %execute runs in parallel
    for i = 1:nruns
        args{1,2}.lookup = args{1,2}.lookup(:,vperm(i,:)); %Permute
        thiscall=args; %need to make a local version of args for parpool to work
        thiscall{3} = seeds(i,:); %third argument is seed
        costs(i) = ClockWalking_3D_discountedEv_stats(thiscall{:});
    end
    
    totcost=sum(costs);
end