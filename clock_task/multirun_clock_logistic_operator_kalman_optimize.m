function [totcost, costs, seeds] = multirun_clock_logistic_operator_kalman_optimize(nruns, runseed, args, vperm)
%use multiple runs of data to identify optimal parameters for logistic
%operator.
rng(runseed);
seeds=randi([1 500], nruns, 4);
costs=zeros(nruns, 1);



%execute runs in parallel --unfortunately I had to change this 
%since I'm guessing matlab can't split assignment operations in parallel

for i = 1:nruns
    args{1,9}.lookup = args{1,9}.lookup(:,vperm(i,:)); %Permute
    thiscall=args; %need to make a local version of args for parpool to work
    thiscall{2} = seeds(i,:); %second argument is seed
    costs(i) = clock_logistic_operator_kalman_optimize(thiscall{:});
end

totcost=sum(costs);
end