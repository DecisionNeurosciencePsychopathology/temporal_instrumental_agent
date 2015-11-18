%Script to optimize parameters of RL models and calculate costs over a
%specified number of runs

%%
%basic structure
%cd '/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/optimality_testing';
agents = initialize_agents_struct;
%agentnames = fieldnames(agents);
%nagents = length(agentnames);
nagents = length(agents);
noptim = 40; %optimize parameters 40x
ncosts = 1000;
ntrials = 100;
nbasis = 24;
ntimesteps = 500;

addpath('../');

%two current approaches: either 'typical' (IEV, DEV, QUADUP, etc.) or 'sinusoid' (phase-shifted random sinusoid)
whichopt=getenv('whichopt');
if strcmpi(whichopt, ''), whichopt='sinusoid'; end

if strcmpi(whichopt, 'typical')
    fprintf('Optimizing parameters using the typical contingencies (IEV, DEV, etc., plus ALL)\n');
    load('mastersamp.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
    condnames=fieldnames(mastersamp); %by default, optimize over all contingencies in lookup
    ncond=length(condnames); %how many conditions/contingencies
    nruns=10; %10 runs per contingency

    %truncate mastersamp to the number of trials used here to save RAM
    for i = 1:length(condnames)
        mastersamp.(condnames{i}).lookup = mastersamp.(condnames{i}).lookup(:,1:ntrials);
    end

    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    optmat = cell(length(condnames),1);
    for i = 1:ncond
        clear row; %need variable to be deleted for array of structs below to work
        for j = 1:nruns
            tmp = mastersamp.(condnames{i});
            tmp.name = condnames{i}; %for keeping track of what contingency this is (esp. when they mix in ALL)
            tmp.lookup = tmp.lookup(1:ntimesteps, randperm(size(tmp.lookup,2))); %randomly permute columns
            tmp.ev = tmp.ev(1:ntimesteps);
            tmp.sample = tmp.sample(1:ntimesteps);
            %tmp.perm = randperm(size(tmp.lookup,2)); %just for checking that this works
            %tmp.lookup = tmp.lookup(:, tmp.perm); %randomly permute columns
            row(j) = tmp; %add to a vector of structs
        end
        optmat{i} = row;
    end

    %concatenate optmat to obtain an ALL condition
    allcond = horzcat(optmat{:});
    optmat{length(optmat) + 1} = allcond;
    condnames{length(condnames) + 1} = 'ALL';
    ncond = ncond + 1; %add all condition

    %optmat is now a cell vector where each element is an vector of structs for the lookups to be used for each run.
    %the final element contains all of the preceding conditions. For the five contingencies and 10 runs, this amounts
    %to a 58MB object that would need to be passed along to parpool workers.
elseif strcmpi(whichopt, 'sinusoid')
    fprintf('Optimizing parameters using a set of phase-shifted random sinusoids\n');
    
    %New tack: load optmat from random sinusoid function treated as one "condition"
    %Similar to the 'ALL' optimization, but with contingencies being semi-random
    load('sinusoid_optmat.mat');
    ncond=1; %optmat only has one element
    
    %truncate to the number of trials for optimization
    for i = 1:length(optmat{1})
        optmat{1}(i).lookup = optmat{1}(i).lookup(1:ntimesteps,1:ntrials);
        optmat{1}(i).ev = optmat{1}(i).ev(1:ntimesteps);
        optmat{1}(i).sample = optmat{1}(i).sample(1:ntimesteps);
    end
elseif strcmpi(whichopt, 'sinusoidsingle')
    %test whether a single permuted sinusoid duplicated many times is stable in optimization
    ncond=1; %optmat has only one element
    nruns=20; %20 runs of a single contingency

    fprintf('Optimizing parameters using a single sinusoid variant permuted %d times\n', nruns);
    load('sinusoid_optmat.mat');
    
    %grab the first element for permutation (NB: this looks like a negative quadratic)
    sinmaster = optmat{1}(1);
    sinmaster.lookup = sinmaster.lookup(1:ntimesteps,1:ntrials); %truncate to the number of trials for optimization
    sinmaster.ev = sinmaster.ev(1:ntimesteps);
    sinmaster.sample = sinmaster.sample(1:ntimesteps);
    
    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    optmat = cell(1,1);
    
    clear row; %need variable to be deleted for array of structs below to work
    for j = 1:nruns
        tmp = sinmaster;
        tmp.lookup = tmp.lookup(:, randperm(size(tmp.lookup,2))); %randomly permute columns
        %tmp.perm = randperm(size(tmp.lookup,2)); %just for checking that this works
        %tmp.lookup = tmp.lookup(:, tmp.perm); %randomly permute columns
        row(j) = tmp; %add to a vector of structs
    end
    optmat{1} = row;
end

%setup random number seeds for each run of data based on a consistent starting seed
rng(888); %an arbitrary consistent starting point for setting seeds that control various probabilistic processes within the agent script
for i = 1:length(optmat)
    for j = 1:length(optmat{i})
        optmat{i}(j).seeds = randi([1 500], 1, 5); %5 random seeds between 1 and 500 per run of data
    end
end

%copy shared optimization parameters to each agent
for i = 1:nagents
    agents(i).nbasis = nbasis;
    agents(i).ntimesteps = ntimesteps;
    agents(i).ntrials = ntrials;
    %agents(i).runseed = 888; %seeds rngs inside runs (deprecated to speed up and maintain greater control)
    
    %agents(i).opt_params = NaN(noptim, ncond, agents(i).k);
    %agents(i).opt_costs = NaN(noptim, ncond);
    %agents(i).opt_conds = condnames;
    
    %agents.(agentnames{i}).opt_params = NaN(noptim, ncond, agents.(agentnames{i}).k);
    %agents.(agentnames{i}).opt_costs = NaN(noptim, ncond);
    %agents.(agentnames{i}).opt_conds = condnames;
end

clear mastersamp sinmaster tmp row allcond i j;

%so, MATLAB requires that for arrays defined outside of loop, but operated on within the parfor,
%the first-level index must be a variant of the parfor counter (i).
%this essentially kills the possibility of loop unrolling because i is never used directly as a first-level index!

ncpus=getenv('matlab_cpus');
if strcmpi(ncpus, '')
    ncpus=40;
    fprintf('defaulting to 40 cpus because matlab_cpus not set\n');
else
    ncpus=str2double(ncpus);
end

%to speed up optimization, have multiple matlab instances that each handle processing for only one agent
whichagent=getenv('which_agent');
if ~strcmpi(whichagent, '')
    pool=['local', whichagent]; %each instance needs its own parpool
    whichagent=str2double(whichagent);
    nagents=1;
    agents=agents(whichagent); %subset only this agent
    system(['touch output/', agents.name, '_running']);
else
    pool='local';
end

%so, it seems that there are two issues here: 1) does the optimizer obtain the same parameters for the same input
% and 2) are optimal parameters stable for slightly different reinforcement histories?
% The first question is how this is setup. The second would require permuting the optmat within each replication.

costs=NaN(noptim, nagents, ncond);
pars=cell(noptim, nagents, ncond);
poolobj=parpool(pool,ncpus);
parfor i = 1:noptim
%for i = 1:noptim

    ipars=cell(nagents, ncond);
    icosts=NaN(nagents, ncond);
    
    for a = 1:nagents
        for c = 1:ncond
            [ipars{a, c}, icosts(a, c)] = GAoptimize(agents(a), optmat{c});
        end
    end
   
    costs(i, :, :) = icosts;
    pars(i, :, :) = ipars;
    
    %save interim progress
    if nagents==1
        parsave(['output/optim_', whichopt, '_', num2str(i), '_', agents(1).name, '.mat'], icosts, ipars, agents); %single agent output
    else
        parsave(['output/optim_', whichopt, '_', num2str(i), '_all.mat'], icosts, ipars, agents); %all agents output
    end
end
delete(poolobj);

if nagents==1
    save(['output/optimize_output_', whichopt, '_', agents(1).name, '.mat'], 'costs', 'pars', 'agents', 'optmat');
else
    %all agents output
    save('output/optimize_output_', whichopt, '_all.mat', 'costs', 'pars', 'agents', 'optmat');
end

%%beautiful, but impossible in matlab. leaving here to cry over

%optimize parameters in parallel.
%need to unroll loop over optimization runs and agents for maximal speed.
%matlab does not allow nested parfor loops, so we need parfor to span out to noptim x nagents.
% nloops = noptim*nagents;
% 
% parfor i = 1:nloops
%     %divide into blocks of length noptim (100)
%     if mod(i-1, noptim) == 0
%         a = (i-1)/noptim + 1;
%     end
%     
%     optnum = mod(i-1, noptim) + 1; %1..100
%     %optresults = NaN(ncond,agents(a).k);
%     optcosts = NaN(ncond,1);
%     
%     %loop over contingencies, optimizing parameters
%     for j = 1:ncond
%         [optresults(j,:), optcosts(j)] = GAoptimize(agents(a), optmat{j});
%     end
%     
%     agents(a).opt_params(optnum,:,:) = optresults; %assign opimal parameters for all conditions of the current optimization to struct
%     agents(a).opt_costs(optnum,:) = optcosts; 
% end
