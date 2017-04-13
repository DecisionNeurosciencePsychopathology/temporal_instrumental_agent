%Script to optimize parameters of clock RL models over a specified number of runs

agents = initialize_agents_struct;
%agentnames = fieldnames(agents);
%nagents = length(agentnames);
nagents = length(agents);
ntrials = 60;
nbasis = 24;
ntimesteps = 500;

addpath('../');

noptim = getenv('noptim');
if strcmpi(noptim, ''), noptim=100;
else noptim=str2double(noptim); end

%two current approaches: either 'typical' (IEV, DEV, QUADUP, etc.) or 'sinusoid' (phase-shifted random sinusoid)
whichopt=getenv('whichopt');
if strcmpi(whichopt, ''), whichopt='sinusoid'; end

%whether to debug
debug=getenv('debug');
if strcmpi(debug, ''), debug=0;
else debug=str2double(debug); end

%whether to fix the softmax beta during optimization (and if so, what value)
fixbeta=getenv('fixbeta');
if strcmpi(fixbeta, ''), fixbeta=0;
else fixbeta=str2double(fixbeta); end

%whether to fix the prop_spread beta during optimization (and if so, what value)
fixps=getenv('fixps');
if strcmpi(fixps, ''), fixps=0;
else fixps=str2double(fixps); end

%number of permuted runs per condition to test (passed to multirun optimizer)
runspercond=getenv('runspercond');
if strcmpi(runspercond, ''), runspercond=10; %default to 10 runs per condition
else runspercond=str2double(runspercond); end

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

rng(621); %an arbitrary consistent starting point for random permutation of contingencies below
if strcmpi(whichopt, 'typical')
    fprintf('Optimizing parameters using the typical contingencies (IEV, DEV, etc., plus ALL)\n');
    load('mastersamp.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
    condnames=fieldnames(mastersamp); %by default, optimize over all contingencies in lookup
    ncond=length(condnames); %how many conditions/contingencies

    %truncate mastersamp to the number of trials used here to save RAM
    for i = 1:length(condnames)
        mastersamp.(condnames{i}).lookup = mastersamp.(condnames{i}).lookup(:,1:ntrials);
    end

    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    optmat = cell(length(condnames),1);
    for i = 1:ncond
        clear row; %need variable to be deleted for array of structs below to work
        for j = 1:runspercond
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
elseif strcmpi(whichopt, 'allequate')
    %Aug2016: four core contingencies (IEV, DEV, CEV, CEVR) with equal EV AUCs (to weigh equally in optimization)
    fprintf('Optimizing parameters for all four Frank contingencies with equal EV AUC\n');
    load('mastersamp_equateauc.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
    condnames={'IEV', 'DEV', 'CEV', 'CEVR'}; %by default, optimize over all contingencies in lookup
    ncond=length(condnames); %how many conditions/contingencies

    %truncate mastersamp_equateauc to the number of trials used here to save RAM
    for i = 1:length(condnames)
        mastersamp_equateauc.(condnames{i}).lookup = mastersamp_equateauc.(condnames{i}).lookup(:,1:ntrials);
    end

    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    optmat = cell(length(condnames),1);
    for i = 1:ncond
        clear row; %need variable to be deleted for array of structs below to work
        for j = 1:runspercond
            tmp = mastersamp_equateauc.(condnames{i});
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
    %optmat{length(optmat) + 1} = allcond;
    %condnames{length(condnames) + 1} = 'ALL';
    %ncond = ncond + 1; %add all condition

    optmat = {allcond};
    condnames = {'ALL'};
    ncond = 1;
    
elseif ismember(whichopt, {'IEV', 'DEV', 'IEVLINPROB', 'DEVLINPROB', 'QUADUP', 'QUADUPOLD', 'QUADDOWN'})
    fprintf('Optimizing parameters using a single contingency %s with %d runs\n', whichopt, runspercond);
    load('mastersamp.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
    condnames=fieldnames(mastersamp); %by default, optimize over all contingencies in lookup
    ncond=1; %only a single condition being optimized

    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    optmat = cell(1,1);
    mastersamp.(whichopt).lookup = mastersamp.(whichopt).lookup(:,1:ntrials); %truncate to number of trials
    
    clear row; %need variable to be deleted for array of structs below to work
    for j = 1:runspercond
        tmp = mastersamp.(whichopt);
        tmp.name = whichopt; %for keeping track of what contingency this is (esp. when they mix in ALL)
        tmp.lookup = tmp.lookup(1:ntimesteps, randperm(size(tmp.lookup,2))); %randomly permute columns
        tmp.ev = tmp.ev(1:ntimesteps);
        tmp.sample = tmp.sample(1:ntimesteps);
        %tmp.perm = randperm(size(tmp.lookup,2)); %just for checking that this works
        %tmp.lookup = tmp.lookup(:, tmp.perm); %randomly permute columns
        row(j) = tmp; %add to a vector of structs
    end
    optmat{1} = row;
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

elseif strcmpi(whichopt, 'doublebump')
    fprintf('Optimizing parameters using a set of two-bump contingencies (low and high)\n');
    
    %load optmat that contains a one-element cell array with 60 varying double bumps
    load('doublebump_optmat.mat');
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
    
    fprintf('Optimizing parameters using a single sinusoid variant permuted %d times\n', runspercond);
    load('sinusoid_optmat.mat');
    
    %grab the first element for permutation (NB: this looks roughly like a negative quadratic)
    sinmaster = optmat{1}(1);
    sinmaster.lookup = sinmaster.lookup(1:ntimesteps,1:ntrials); %truncate to the number of trials for optimization
    sinmaster.ev = sinmaster.ev(1:ntimesteps);
    sinmaster.sample = sinmaster.sample(1:ntimesteps);
    sinmaster.prb = sinmaster.prb(1:ntimesteps);
    sinmaster.mag = sinmaster.mag(1:ntimesteps);
    
    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    optmat = cell(1,1);
    
    clear row; %need variable to be deleted for array of structs below to work
    for j = 1:runspercond
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
        optmat{i}(j).firstrt = randi([1 ntimesteps], 1, 1); %set random first RT for choice rule.
    end
end

%copy shared optimization parameters to each agent
for i = 1:nagents
    agents(i).nbasis = nbasis;
    agents(i).ntimesteps = ntimesteps;
    agents(i).ntrials = ntrials;
    agents(i).debug = debug; %debug optimization (don't actually run GA)
    
    %handle fixed beta by removing it from agent structure when beta is a relevant parameter
    if fixbeta > 0
        if any(ismember(agents(i).parnames, 'beta'))
            bpos = ismember(agents(i).parnames, 'beta');
            agents(i).init_params = agents(i).init_params(~bpos);
            agents(i).lower_bounds = agents(i).lower_bounds(~bpos);
            agents(i).upper_bounds = agents(i).upper_bounds(~bpos);
            agents(i).k = agents(i).k - 1;
            agents(i).parnames = agents(i).parnames(~bpos);
            agents(i).fixbeta = fixbeta;
        end
    else
        agents(i).fixbeta=0;
    end
    
    %handle fixed prop_spread by removing it from agent structure
    if fixps > 0
        if any(ismember(agents(i).parnames, 'prop_spread'))
            pspos = ismember(agents(i).parnames, 'prop_spread');
            agents(i).init_params = agents(i).init_params(~pspos);
            agents(i).lower_bounds = agents(i).lower_bounds(~pspos);
            agents(i).upper_bounds = agents(i).upper_bounds(~pspos);
            agents(i).k = agents(i).k - 1;
            agents(i).parnames = agents(i).parnames(~pspos);
            agents(i).fixps = fixps;
        end
    else
        agents(i).fixps=0;
    end
end

clear mastersamp mastersamp_equateauc sinmaster tmp row allcond i j;

%so, MATLAB requires that for arrays defined outside of loop, but operated on within the parfor,
%the first-level index must be a variant of the parfor counter (i).
%this essentially kills the possibility of loop unrolling because i is never used directly as a first-level index!

%so, it seems that there are two issues here: 1) does the optimizer obtain the same parameters for the same input
% and 2) are optimal parameters stable for slightly different reinforcement histories?
% The first question is how this is setup. The second would require permuting the optmat within each replication.

if fixbeta > 0, betastr=['_beta' num2str(fixbeta)];
else betastr=''; end

if fixps > 0, psstr=['_ps' num2str(fixps)];
else psstr=''; end

costs=NaN(noptim, nagents, ncond);
pars=cell(noptim, nagents, ncond);
%poolobj=parpool(pool,ncpus);
poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
ifiles=cell(1, noptim);

try
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
            ifiles{i} = ['output/optim_', whichopt, betastr, psstr, '_', num2str(i), '_', agents(1).name, '.mat'];
        else
            ifiles{i} = ['output/optim_', whichopt, betastr, psstr, '_', num2str(i), '_all.mat']; %all agents output
        end
        
        parsave(ifiles{i}, icosts, ipars, agents); %single agent output
    end
catch err
    disp('error in optimization. killing parpool');
    delete(poolobj);
    rethrow(err);
end

delete(poolobj);

if nagents==1
    save(['output/optimize_output_', whichopt, betastr, psstr, '_', agents(1).name, '.mat'], 'costs', 'pars', 'agents', 'optmat');
else
    %all agents output
    save(['output/optimize_output_', whichopt, betastr, psstr, '_all.mat'], 'costs', 'pars', 'agents', 'optmat');
end

%if we have made it here, then all optimizations have completed successfully, so cleanup the interim files
ifiles = ifiles(~cellfun('isempty', ifiles)); %remove any empty elements
delete(ifiles{:}) %delete all intermediate files

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
