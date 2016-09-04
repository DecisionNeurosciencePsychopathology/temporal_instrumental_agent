%% Vet TC identifiability challenges against Michael Frank's TC Feedback Aug2016
addpath('../'); %for the RewFunction code
addpath('../optimality_testing'); %for the TC forward code

%Simulation parameters from Michael Frank email 19Aug2016
% K = 1500;
% lambda = 0.2;
% epsilon = 3000;
% alphaG = 1;
% alphaN = .5;
% rho = 2500; % (or 0) 
% nu = .3;

format shortG; %easier to look at

%params = [lambda, epsilon, alphaG, alphaN, K, nu, rho]; %in order from original TC (used in forward)
parbounds = [ ...
    0, 1; ... %lambda
    0, 20000; ... %epsilon
    0, 5; ... %alphaG
    0, 5; ... %alphaN
    0, 5000; ... %K    
    0, 10; ... %nu -- unclear how to set this range
    0, 10000; ... %rho
    ];

%fixed K, lambda, nu (1500, 0.2, 0.2 per MJF identifiability email)
% parbounds = [ ...
%     0.2, 0.2; ... %lambda
%     0, 50000; ... %epsilon
%     0.1, 3; ... %alphaG
%     0.1, 3; ... %alphaN
%     1500, 1500; ... %K  
%     0.2, 0.2; ... %nu -- unclear how to set this range
%     0, 10000; ... %rho
%     ];

nreps = 1000;
npars = size(parbounds, 1);

rng(101); %fixed parsets
parsets = NaN(nreps, npars);

for i = 1:nreps
    for j = 1:npars
        parsets(i,j) = parbounds(j, 1) + (parbounds(j,2) - parbounds(j, 1)).*rand(1,1);
    end
end

%relax bounds on RT here to rule out possibility that nonidentiability is due to truncated RT range
rtbounds = [-30000 30000]; %define valid choice RTs (since agent is choosing freely)

%Per discussion with MJF, allow truly random draws from distribution over simulations.
%This eliminates perfect repeatability across simulations, but since we are averaging and smoothing
%across 'subjects,' the results should be convergent.
userngseed = 0;
ntrials = 50;

%% Simulate data (forward model)
%complete one run of each contingency (order permuted within replications loop)
conds = {'DEV', 'IEV', 'CEV', 'CEVR'};

%allRTpred = NaN(nreps, length(conds), ntrials);
%allRewpred = NaN(nreps, length(conds), ntrials);

clear simstruct

for i = 1:nreps
   condperm = conds(randperm(length(conds)));
   %setup priors for initial run
   priors.V = 0;
   priors.Go = 0;
   priors.NoGo = 0;
   %priors.AvgRT = 2000; %average RT
   %priors.FirstRT = 2000; %first choice

   %Gaussian variation in mean RT
   priors.AvgRT = normrnd(2000, 400, 1);
   
   %Gaussian variation in first choice
   priors.FirstRT = normrnd(2000, 400, 1);
      
   for j = 1:length(condperm)
       [~, RTpred, ret] = TC_Alg_forward(parsets(i,:), priors, condperm{j}, -1, ntrials, rtbounds); %-1 for rngseed indicates not to use repeatable rng
       simstruct(i,j).rt = RTpred;
       simstruct(i,j).rew = ret.rew;
       simstruct(i,j).pars = parsets(i,:);
       simstruct(i,j).cond = condperm{j};
       simstruct(i,j).priors = priors;

       %allRTpred(i, ismember(conds, condperm{j}), :) = RTpred; %find undoes the permutation so that indices follow cond order above
       %allRewpred(i, find(ismember(conds, condperm{j})), :) = ret.rew;
       %RTsmooth = smooth(RTpred, span, smoother);
       %allRTsmooth(i, find(ismember(conds, condperm{j})), :) = RTsmooth;
       priors.V = ret.V(length(ret.V)); %carry forward value to next run
   end
end

%average over subjects
% allRTavg = squeeze(mean(allRTpred, 1));
% allRTsmoothavg = squeeze(mean(allRTsmooth, 1));

%plot results
% figure(1); plot(allRTavg', 'LineWidth',2);
% legend(conds); xlabel('Trial'); ylabel('Average RT across replications'); title(sprintf('TC results with rho = %.0f', rho));

%figure(2); plot(allRTsmoothavg', 'LineWidth',2);
%legend(conds); xlabel('Trial'); ylabel('Average RT across replications'); title(sprintf('Smoothed TC results with rho = %.0f', rho));
%print('TC_sims_withoutrho','-dpng','-r200')

%% Fit simulated data -- recover parameters (identifiability) 

num_start_pts = 8; % number of initial starting points for parameter estimation

%optimizer settings
opts = optimset('fmincon');
opts.LargeScale = 'off';
opts.Algorithm = 'active-set';
opts.Display = 'none';

model = 'noemo'; %vanilla model under old vetted fitting code

init_params = [0.3; 2000; 0.2; 0.2; 1000; 0.1; 0.5; 300]; %lambda, epsilon, alphaG, alphaN, K, nu, expalt (UNUSED), rho
lower_limits = parbounds(:,1);
upper_limits = parbounds(:,2);

%need to add in exp_alt dummy since this is not part of forward algorithm, but is expected in fitting (unused)
lower_limits = [lower_limits(1:6); 0; lower_limits(7)];
upper_limits = [upper_limits(1:6); 0; upper_limits(7)];

fittedpars = NaN(nreps, npars);

ncpus=getenv('matlab_cpus');
if strcmpi(ncpus, '')
    ncpus=40;
    fprintf('defaulting to 40 cpus because matlab_cpus not set\n');
else
    ncpus=str2double(ncpus);
end

poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)

try

%loop over replications, fitting parameters using fmincon with multiple starting points
parfor i = 1:nreps
    %pass in a struct array to TC_minSE containing one element for each block
        
    %core fitting function -- returns results of length num_start_pts (8) above.
    %These contain fit estimates for each starting point.
    %Then identify the best-fitting output for use in analyses.
    [params, SE, exitflag, xstart] = rmsearch(@(params) TC_minSE(params, simstruct(i,:), model), 'fmincon', init_params, ...
        lower_limits, upper_limits, 'initialsample', num_start_pts, 'options', opts);
    
    bestpars = params(find(SE == min(SE), 1 ),:);
    %re-run TC_alg for all blocks with optimal parameters
    [totalSqErr, ret_all] = TC_minSE(bestpars, simstruct(i,:), model);
    
    fittedpars(i, :) = bestpars([1:6, 8]);
      
end

catch err
    disp('error in optimization. killing parpool');
    delete(poolobj);
    rethrow(err);
end

delete(poolobj);

save('tc_identifiability.mat', 'parsets', 'fittedpars', 'parbounds', 'rtbounds', 'lower_limits', 'upper_limits', 'init_params', 'simstruct');
