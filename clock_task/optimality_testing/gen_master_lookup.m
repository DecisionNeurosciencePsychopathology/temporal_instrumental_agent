%This script creates a master lookup table for each condition (IEV, DEV, QUADUP, IEVLINPROB, DEVLINPROB).
%The rows represent time steps (10ms bins) and the columns represent successive draws from the contingency.
%NOTE: the TD agents sample at 100ms intervals I believe so when
%sampling occurs be sure to reduce (or remap) the master table so it can
%sample properly. Is this inherently a hinderence on the Q models???

%The idea is to have the same series of outcomes for an agent who samples the same RTs.
%So if SKEPTIC chooses 3500ms 4 times and Q-learning chooses 3500ms 4 times, both receive the same series
%of outcomes. This eliminates this source of variability in the costs.

addpath('../');

ntimesteps = 500; %in 10ms bins (1-5s)
conds = {'IEV' 'DEV' 'QUADUP' 'QUADUPOLD' 'QUADDOWN', 'IEVLINPROB' 'DEVLINPROB'};
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
    mastersamp.(conds{i}).sample = zeros(1,ntimesteps); %keeps track of how many times a timestep has been sampled by agent
    mastersamp.(conds{i}).ev = zeros(1,ntimesteps);
    
    for j = 1:ntimesteps
        [~, mastersamp.(conds{i}).ev(j), mastersamp.(conds{i}).prb(j), mastersamp.(conds{i}).mag(j)] = RewFunction(j*10, conds{i}, 0, 5000);
        for k = 1:ntrials
            [mastersamp.(conds{i}).lookup(j,k)] = RewFunction(j*10, conds{i}, 0, 5000);
        end
    end
    
    %parsave(['mastersamp_' conds{i} '.mat'], mastersamp);
end
save('mastersamp.mat', 'mastersamp');

%New approach: generate variants of time-varying contingencies that do not prefer sampling at the edge.
%Use sinusoidal functions to generate a continuous contingency that can be shifted in time to maintain a constant EV AUC (identical costs)
%while shifting the optimal RT. Then optimize a subset of these possible contingencies (sampled at random)

ntimesteps=500;

ev = 10*sin(2*pi*(1:ntimesteps).*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps);
ev = ev + abs(min(ev)) + 10;
prb = 25*cos(2*pi*(1:ntimesteps).*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps);
prb_max=0.7;
prb_min=0.3;
prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min;

%simpler version without substantial high-frequency oscillation
% ev = 10*sin(2*pi*(1:ntimesteps).*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps);
% ev = ev + abs(min(ev)) + 10;
% prb = 25*cos(2*pi*(300-(1:ntimesteps)).*1/ntimesteps); %+ 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps);
% prb_max=0.7;
% prb_min=0.3;
% prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min;

%mag = mag + abs(min(mag)) + 10;
%prb = evi./magi;
%plot(1:ntimesteps, mot1, type="l")

%figure(1); clf;
allshift = NaN(ntimesteps, ntimesteps, 3);
conds = 1:ntimesteps;
for i = 1:ntimesteps
  shift=[i:ntimesteps 1:(i-1)];
  evi = ev(shift);
  prbi = prb(shift);
  
  allshift(i,:,1) = evi;
  allshift(i,:,2) = prbi;
  allshift(i,:,3) = evi./prbi;
  
  %subplot(3,1,1); plot(1:ntimesteps, evi)
  %subplot(3,1,2); plot(1:ntimesteps, prbi);
  %subplot(3,1,3); plot(1:ntimesteps, evi./prbi);
  %pause(0.02)
  
end

rng(102); %fix seed for pulling reward probabilities
ntrials = 500; %maximum number of trials that could be used for testing this contingency

%randomly sample 60 of the possible 500 contingencies without replacement
keep = randsample(1:ntimesteps, 60);

%for consistency with prior simulations, keep sinusoid 366 as first variant (this is roughly quad down with max at 250 and prb max at 150)
keep(1) = 366;

optmat=cell(1,1);
for k = 1:length(keep)
    thisCont=[];
    thisCont.name = ['sinusoid' num2str(keep(k))];
    thisCont.sample = zeros(1, ntimesteps); %keeps track of how many times a timestep has been sampled by agent
    thisCont.lookup = zeros(ntimesteps, ntrials); %lookup table of timesteps and outcomes
    thisCont.ev = allshift(keep(k),:,1);
    thisCont.prb = allshift(keep(k),:,2);
    thisCont.mag = allshift(keep(k),:,3);
    
    rvec = rand(ntimesteps, ntrials);
    for t = 1:ntrials
        thisCont.lookup(:,t) = (allshift(keep(k),:,2) > rvec(:,t)') .* allshift(keep(k),:,3);
    end
    
    optmat{1}(k) = thisCont;
end

save('sinusoid_optmat.mat', 'optmat');
