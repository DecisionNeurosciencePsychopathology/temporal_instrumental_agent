function [totcost, costs, seeds] = multirun_skeptic_fitsubsubect_all_models(nruns, runseed, args)
%use multiple runs of data to identify optimal parameters for logistic
%operator.
rng(runseed);
seeds=randi([1 500], nruns, 4);
costs=zeros(nruns, 1);

%default range seeds for ref
%rngseeds=[98 83 66 10];

%execute runs in parallel --unfortunately I had to change this 
%since I'm guessing matlab can't split assignment operations in parallel

%Try parfor if it breaks go back to for, remove automatic parfor loop in
%preferences, bottom left cornor in main GUI -> prefs
for i = 1:nruns
    thiscall=args; %need to make a local version of args for parpool to work
    thiscall{4} = seeds(i,:); %second argument is seed
    costs(i) = skeptic_fitsubject_all_models(thiscall{:});
end

totcost=sum(costs);
end