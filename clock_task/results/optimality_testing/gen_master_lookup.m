%This script creates a master lookup table for each condition (IEV, DEV, QUADUP, IEVLINPROB, DEVLINPROB).
%The rows represent time steps (10ms bins) and the columns represent successive draws from the contingency.
%NOTE: the TD agents sample at 100ms intervals I believe so when
%sampling occurs be sure to reduce (or remap) the master table so it can
%sample properly. Is this inherently a hinderence on the Q models???

%The idea is to have the same series of outcomes for an agent who samples the same RTs.
%So if SKEPTIC chooses 3500ms 4 times and Q-learning chooses 3500ms 4 times, both receive the same series
%of outcomes. This eliminates this source of variability in the costs.

ntimesteps = 500; %in 10ms bins (1-5s)
conds = {'IEV' 'DEV' 'QUADUP' 'IEVLINPROB' 'DEVLINPROB'};
ntrials = 500; %number of draws for each RT where each draw becomes a column in the lookup matrix

%set up seeds
%global rew_rng_state;
%rew_rng_seed=71; %seed used to populate outcomes 
%rng(rew_rng_seed);
%rew_rng_state=rng;

mastersamp=[];
for i = 1:length(conds)
    %Initalize structure
    mastersamp.(conds{i}).lookup = zeros(ntimesteps,ntrials); %lookup table of timesteps and outcomes
    mastersamp.(conds{i}).sample = ones(1,ntrials); %keeps track of how many times a timestep has been sampled by agent
    mastersamp.(conds{i}).ev = zeros(1,ntrials);
    
    for j = 1:ntimesteps
        [~, mastersamp.(conds{i}).ev(j)] = RewFunction(j*10, conds{i}, 0, 5000);
        for k = 1:ntrials
            [mastersamp.(conds{i}).lookup(j,k)] = RewFunction(j*10, conds{i}, 0, 5000);
        end
    end
    
    %parsave(['mastersamp_' conds{i} '.mat'], mastersamp);
end
save('mastersamp.mat', 'mastersamp');

