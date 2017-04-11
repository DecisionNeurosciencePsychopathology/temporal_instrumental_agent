%plan of attack
%1) generate 20 datasets at each of 5 optimal parameter sets for each model (n = 100). Trials = 50
%2) fit each of the 100 datasets using all other models (using rmsearch or VBA?)
%3) compute relative fit statistics for data-generating model versus fitting model

%cd '/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing/output';
%cd '/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/optimality_testing';
%addpath('../');
basedir='/Users/michael/Data_Analysis/temporal_instrumental_agent';
addpath([basedir, '/clock_task']);
addpath([basedir, '/clock_task/optimality_testing']);
addpath([basedir, '/clock_task/vba']);
optfiles = glob([basedir, '/clock_task/optimality_testing/output/optimize_output_sinusoid*.mat']);
%optfiles = glob('output/optimize_output_allequate*.mat');
%optfiles = glob('output/optimize_output_sinusoid_fixedLR_decay*.mat');

optfiles=optfiles([1 4 5 9 11]);

nreps = 40;
nbest = 3; %number of parameter sets to test
%nbest=1; %just for decay tests
ntimesteps=500;
%runlengths = [30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
%runlengths = [30 40 60 110];
runlengths = 100; %just one run for confusion, but long enough to see convergence
cliffpullback=30;
plots=false;
%target='allequate';
target='sinusoid'; %contingency for simulation

if strcmpi(target, 'sinusoid')
    %draw a random sample from the phase-shifted sinusoids for each replication
    ev = 10*sin(2*pi*(1:ntimesteps).*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps);
    ev = ev + abs(min(ev)) + 10;
    prb = 25*cos(2*pi*(1:ntimesteps).*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps);
    prb_max=0.7;
    prb_min=0.3;
    prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min;
    
    allshift = NaN(ntimesteps, ntimesteps, 3);
    conds = 1:ntimesteps;
    for i = 1:ntimesteps
        shift=[i:ntimesteps 1:(i-1)];
        evi = ev(shift);
        prbi = prb(shift);
        
        allshift(i,:,1) = evi;
        allshift(i,:,2) = prbi;
        allshift(i,:,3) = evi./prbi;
    end
    
    %randomly sample with replacement the desired number of contingencies
    keep = randsample(1:ntimesteps, nreps, true);
    
    rng(625); %fix seed for pulling reward probabilities
    
    clear optmat;
    for k = 1:length(keep)
        thisCont=[];
        thisCont.name = ['sinusoid' num2str(keep(k))];
        thisCont.sample = zeros(1, ntimesteps); %keeps track of how many times a timestep has been sampled by agent
        thisCont.lookup = zeros(ntimesteps, max(runlengths)); %lookup table of timesteps and outcomes
        thisCont.ev = allshift(keep(k),:,1);
        thisCont.prb = allshift(keep(k),:,2);
        thisCont.mag = allshift(keep(k),:,3);
        
        rvec = rand(ntimesteps, max(runlengths));
        for t = 1:max(runlengths)
            thisCont.lookup(:,t) = (allshift(keep(k),:,2) > rvec(:,t)') .* allshift(keep(k),:,3);
        end
        
        optmat(k) = thisCont;
    end
    
elseif strcmpi(target, 'allequate')
    %Aug2016: four core contingencies (IEV, DEV, CEV, CEVR) with equal EV AUCs (to weigh equally in optimization)
    load('mastersamp_equateauc.mat'); %sampling lookup for all contingencies (maintain identical rewards for identical choices)
    condnames={'IEV', 'DEV', 'CEV', 'CEVR'}; %by default, optimize over all contingencies in lookup
    ncond=length(condnames); %how many conditions/contingencies
    ntimesteps = 500;
    
    repspercond = floor(nreps/length(condnames));
    condmat = repmat(condnames, repspercond, 1);
    nrep = cell(1, nreps);
    for i = 1:ncond
        for j = 1:repspercond
            nrep{(i-1)*repspercond + j} = num2str(j);
        end
    end
    %condnames = strcat(condmat(:)', nrep);
    condrep = condmat(:)';
    
    %generate master lookup, permuting the lookup columns for each run to represent different draws from the contingency
    %handle permutation outside of optimization for speed
    clear optmat;
    for i = 1:ncond
        for j = 1:repspercond
            tmp = mastersamp_equateauc.(condnames{i});
            tmp.name = condnames{i}; %for keeping track of what contingency this is (esp. when they mix in ALL)
            tmp.lookup = tmp.lookup(1:ntimesteps, randperm(size(tmp.lookup,2))); %randomly permute columns
            tmp.ev = tmp.ev(1:ntimesteps);
            tmp.sample = tmp.sample(1:ntimesteps);
            optmat((i-1)*repspercond + j) = tmp; %add to a vector of structs
        end
    end
end

nagents=length(optfiles);

priors.V = 0; %don't give agent any insight into previous values (no SCEPTIC receives this boost)
priors.Go = 0;
priors.NoGo = 0;

%fly through the contingencies (looks good)
%for i = 1:length(optmat)
%    pause(0.08);
%    subplot(3,1,1);
%    plot(optmat(i).prb);
%    subplot(3,1,2);
%    plot(optmat(i).mag);
%    subplot(3,1,3);
%    plot(optmat(i).ev);
%end

%storage structure for output
oresults=struct();

%prop_spread const
prop_spread_const = .05;
beta_const = 1.5;


for o = 1:(length(optfiles) + 1)
    if o <= length(optfiles);
        opt = load(optfiles{o});
        a=opt.agents;
        
        %look at top 10% and compare parameters
        %bestqcosts = sort(opt.costs);
        %bestqcosts = bestqcosts(1:nbest);
        %cutoff = quantile(opt.costs, .1); %bottom 10% since lower costs are better
        %cutoff = max(bestqcosts);
        %bestcost_indices = opt.costs <= cutoff;
        %bestcosts = opt.costs(bestcost_indices);
        %[~, ind] = sort(bestcosts);
        %candidates = cell2mat(opt.pars(bestcost_indices));
        %candidates = horzcat(candidates, bestcosts);
        %candidates = sortrows(candidates, size(candidates, 2));
        %candidates = candidates(ind, :); %sort in terms of low to high (good to bad) costs)
        
        
        if strcmpi(a.name, 'fixedLR_decay')
            candidates = [ prop_spread_const beta_const 0.064386 -1.6371; ... %prop_spread, beta, alpha, gamma
                prop_spread_const beta_const 0.052754 -1.1632; ... 
                prop_spread_const beta_const 0.06 -0.5];            
            %candidates = [ 0.035305       1.6889     0.064386      -1.6371; ... %prop_spread, beta, alpha, gamma
            %    0.062697       1.9736     0.052754      -1.1632];
        elseif strcmpi(a.name, 'fixedLR_softmax')
            candidates = [ prop_spread_const beta_const 0.04; ... %prop_spread, beta, alpha
                prop_spread_const beta_const 0.06; ...
                prop_spread_const beta_const 0.08 ];
            %candidates = [      0.04174       1.6862     0.043147; ... %prop_spread, beta, alpha
            %    0.037578      1.9964     0.051624];
        elseif strcmpi(a.name, 'fixed_uv')
            candidates = [ prop_spread_const beta_const .1 .1; ... %prop_spread, beta, alpha, tau
                prop_spread_const beta_const .1 .2; ...
                prop_spread_const beta_const .1 .3 ];            
            %candidates = [ 0.06872       1.4279     0.095586     0.019568; ... %prop_spread, beta, alpha, tau
            %    0.033588      1.1007     0.099512     0.015764];
        elseif strcmpi(a.name, 'kalman_softmax')
            candidates = [ prop_spread_const beta_const; ... %prop_spread, beta
                prop_spread_const beta_const;
                prop_spread_const beta_const ];
            %candidates = [     0.19576       1.7564; ... 
            %    0.22062       1.7318 ];
        elseif strcmpi(a.name, 'kalman_uv_sum_negtau')
            candidates = [ prop_spread_const beta_const .1; ... %prop_spread, beta, tau
                prop_spread_const beta_const .2; ...
                prop_spread_const beta_const .3 ];
            %candidates = [ 0.18117      0.17702     0.057835; ... 
            %    0.1541      0.13649     0.062604];
        end
    else
        a.name = 'null';
        a.nbasis = 24;
        a.ntimesteps = 500;
        a.prop_spread=0; %need this to avoid problem with parameter loading
        a.fixps=1; a.fixbeta=0; %hack: fix prop spread to have the parameter lookup work
        parsq=[];
    end
    
    rng(1040); %consistent starting point for randomization of contingency
    
    %pars = cell2mat(opt.pars);
    %medpars = median(pars);
    %bestpars = pars(find(opt.costs==min(opt.costs)),:);
    
    %v_it=NaN(nagents, nreps, runlengths, ntimesteps);
    %rts = NaN(length(runlengths), nreps, nbest, runlengths);
    %ret = cell(nagents, nreps);
    
    oresults.(a.name).costs = NaN(length(runlengths), nreps, nbest);
    oresults.(a.name).pars  = cell(length(runlengths), nreps, nbest);
    oresults.(a.name).rts   = cell(length(runlengths), nreps, nbest);
    oresults.(a.name).rewards = cell(length(runlengths), nreps, nbest);
    for t = 1:length(runlengths)
        %rts   = NaN(runlengths(t),1);
        
        for p = 1:nreps
            seeds = randi([1 500], 1, 5);
            optmat(p).firstrt = randi([1 500], 1, 1);
            for q = 1:nbest
                %parsq = candidates(q, 1:(size(candidates,2)-1)); %trim last column, which contains cost
                parsq = candidates(q, :); %trim last column, which contains cost
                %oresults.(a.name).pars{t,p,q} = parsq;
                if strcmpi(a.name, 'franktc')
                    %[costs(o,p), rts(o,p,:), ret{o,p}] = TC_Alg_forward(bestpars, priors, optmat(p), seeds, a.ntrials, [0 5000]);
                    [oresults.(a.name).costs(t,p,q), oresults.(a.name).rts{t,p,q}, ~] = TC_Alg_forward(parsq, priors, optmat(p), seeds, runlengths(t), [0 5000]);
                elseif strcmpi(a.name, 'qlearning')
                    a.clock_options.episodeCount = runlengths(t);
                    a.clock_options.ntimesteps = a.ntimesteps/10; %TD models operate on 50 timesteps, not 500, typically
                    [oresults.(a.name).costs(t,p,q), ~, oresults.(a.name).rts{t,p,q}] = ClockWalking_3D_discountedEv_optimize(a.clock_options, optmat(p), seeds, parsq);
                else
                    %[costs(o,p), v_it(o,p,:,:), rts(o,p,:), ret{o,p}] = clock_sceptic_agent(bestpars, a, seeds, optmat(p), ntrials, a.nbasis, a.ntimesteps);
                    %[costs(o,p), ~, rts(o,p,:), ret{o,p}] = clock_sceptic_agent(bestpars, a, seeds, optmat(p), ntrials, a.nbasis, a.ntimesteps);
                    [oresults.(a.name).costs(t,p,q), ~, oresults.(a.name).rts{t,p,q}, ret] = clock_sceptic_agent(parsq, a, seeds, optmat(p), runlengths(t), a.nbasis, a.ntimesteps, cliffpullback);
                    oresults.(a.name).rewards{t,p,q} = ret.rew_i;
                end
            end
        end
    end
end

%save('/Users/michael/TresorSync/confusionsim.mat', 'oresults', 'optmat', 'optfiles'); %'-v7.3'
%this makes the file about 8GB!  'ret' 'v_it',
%save('/Users/michael/bootstrap_costs_at_best10.mat', 'costs', 'rts', 'optmat', 'optfiles'); %'-v7.3'

%fork of sceptic_fit_group_vba


%% choose models to fit
fitmodelnames = {'fixed_decay' 'fixed' 'fixed_uv' 'kalman_softmax' 'kalman_uv_sum'}; %models to fit to each dataset
simmodelnames = fieldnames(oresults);
nsimmodels = length(simmodelnames);
nfitmodels = length(fitmodelnames);
ncpus=40;

%% set parameters
nbasis = 24;
multinomial = 1;
multisession = 0;
fixed_params_across_runs = 1;
fit_propspread = 1;
n_steps = 50;

u_aversion = 1; % allow for uncertainty aversion in UV_sum
saveresults = 1;
graphics = 0; %display fitting

% get ID list
%id = cell(nreps,1);
id = NaN(nreps,1);


%% main loop
L = NaN(nsimmodels,nfitmodels,nbest,nreps);
poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
grp = struct([]);

for i=1:nsimmodels
    simmodel = simmodelnames{i};
    results_dir = ['simdf_scratch/', simmodel];
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    for j=1:nfitmodels
        for k = 1:nbest
            fitmodel = char(fitmodelnames(j));
            nsubjects = size(oresults.(simmodel).rts, 2); %rts is triallength x nreplications x parset, so here 1 x 50 x 2
            
            parfor sub=1:nsubjects
                %write out relevant data to disk
                cname=[];
                [cname{1:runlengths}] = deal(optmat(sub).name);
                df = table(ones(runlengths, 1), (1:runlengths)', cname', (oresults.(simmodel).rts{1, sub, k}.*10)', ... %multiply by 10 to go from 500 -- 5000 scaling
                    oresults.(simmodel).rewards{1, sub, k}', 'VariableNames', {'run', 'trial', 'rewFunc', 'rt', 'score'});
                datafile = sprintf('simdf_scratch/simdf_%s_%d.csv', simmodel,sub);
                writetable(df,datafile);
                %id{sub} = [simmodel, '_', sub];
                id(sub) = sub; %keep id as numeric for now
                                
                fprintf('simmodel: %s, fitmodel: %s, subject: %d \r',simmodel, fitmodel, sub)
                [posterior,out] = clock_sceptic_vba(id(sub), fitmodel,nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, ...
                    n_steps, u_aversion, datafile,saveresults, graphics, results_dir, datafile, optmat(sub)); %pass in the contingency matrix from sim as last param
                
                %             cd(results_dir);
                %             parsave(sprintf('output_%d',id(sub)),posterior,out);
                L(i,j,k,sub) = out.F;
            end
        end
    end
end

delete(poolobj);

save('confusion_out.mat', 'L', 'fitmodelnames', 'simmodelnames', 'oresults', 'optmat');

%cd(group_dir);
%filename = sprintf('SHIFTED_U_grp_L_%d_nbasis%d_nsteps%d_uaversion_not_allModels_fixed_prop_spread_frank_all_subjs',nbasis,n_steps, u_aversion);
%save(filename);