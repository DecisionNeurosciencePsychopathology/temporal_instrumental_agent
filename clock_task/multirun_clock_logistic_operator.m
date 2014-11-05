function [totcost, costs, seeds] = multirun_clock_logistic_operator(nruns, runseed, varargin)
    %use multiple runs of data to identify optimal parameters for logistic
    %operator.
    rng(runseed);
    seeds=randi([1 500], nruns, 2);
    
    
    args=varargin;
    costs=zeros(nruns, 1);
    for i = 1:nruns
        args{2} = seeds(i,:); %second argument is seed
        costs(i) = clock_logistic_operator(args{:});
    end
    
    totcost=sum(costs);
end