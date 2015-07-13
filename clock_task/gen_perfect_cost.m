%Create the perfect cost value from x amount of trials for each condition,
%I believe we can do this via the RewFucntion or by chosing always choosing
%the max value in the pseudo random loop and averaging it. Then you would
%be comparing the agents average to the perfect cost average.

ntimesteps = 500; %in 10ms bins
conds = {'IEV' 'DEV' 'QUADUP'};
ntrials = 200;

%set up seeds--Delete, doesnt matter...
global rew_rng_state;
rew_rng_seed=71;
rng(rew_rng_seed);
rew_rng_state=rng;




%This will calculate the Expected value for the number of time steps
for i = 1:length(conds)
    for k = 1:ntimesteps
        [~,perfect_cost.(conds{i})(k,1)] = RewFunction(k*10,conds{i},1);        
    end
    [perfect_cost.([conds{i} '_pertrial'])] = max(perfect_cost.(conds{i})); %For 1 trial
    [perfect_cost.(conds{i})] = max(perfect_cost.(conds{i}))*ntrials; %For ntrials
end

    

%perfect_cost.reversal = (perfect_cost.IEV/2 + perfect_cost.DEV/2);
%perfect_cost.reversal_pertrial = (perfect_cost.IEV_pertrial + perfect_cost.DEV_pertrial)/2;