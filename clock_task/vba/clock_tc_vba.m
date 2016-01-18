function [posterior,out] = clock_tc_vba(id,multinomial,multisession,fixed_params_across_runs)
%% fits TC model to Clock Task subject data using VBA toolbox
% example call:
% [posterior,out]=clock_sceptic_vba(10638,'fixed',4,1,1,1)
% id:           5-digit subject id in Michael Hallquist's BPD study
% only works with 'fixed' (fixed learning rate SCEPTIC) so far
% n_basis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)

close all

n_t = 400; %400 trials
n_runs = 8;
trialsToFit = 1:n_t;

%% fit as multiple runs
multisession = 1;
% fix parameters across runs
fixed_params_across_runs = 1;

%% u is 2 x ntrials where first row is rt and second row is reward
data = readtable(sprintf('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
rts = data{trialsToFit, 'rt'};
rewards = data{trialsToFit, 'score'};

%% split into conditions/runs
if multisession %improves fits moderately
    options.multisession.split = repmat(n_t/n_runs,1,n_runs);
    %% fix parameters
    if fixed_params_across_runs
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        % allow unique initial values for each run?
        %                  options.multisession.fixed.X0 = 'all';
    end
end

%in Frank TC, n_phi is 5 and n_theta is 2
%number of hidden states (n) is 10
n_states = 10;

dim = struct('n',n_states,'n_theta',2,'n_phi',5, 'n_t', n_t);
y = data{trialsToFit,'rt'}';
priors.a_alpha = Inf;   %infinite precision prior on state noise precision (deterministic system)
priors.b_alpha = 0;
priors.a_sigma = 1;     % Jeffrey's prior (on measurement noise precision)
priors.b_sigma = 1;     % Jeffrey's prior

options.binomial = 0;
options.sources(1) = struct('out',1,'type',0); %treat source/y as Gaussian

%% skip first trial (don't fit t=1)
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;

%fix alphaG and alphaN to 0-1, exp transform, and multiply by 5 inside the function

%treat K as Gaussian with mean in the middle of the interval (2000) and high variance. Try to achieve a truly flat prior (uniform)
%treat lambda as exponential transform 0-1 hack
%treat epsilon as a gamma function (long tail) -- but VBA only allows for Gaussian priors. Use a distributional transform
%treat rho as Gaussian with positive mean and large variance

%for X0, use mean = 0 and sigma = 0 to avoid fitting initial values for hidden states (a_short, b_short, etc.)

%% priors
%theta is alphaG (1), alphaN (2)
priors.muTheta = [-3, -3]; %learning rate prior of .2371: 5/(1+exp(-theta))
priors.SigmaTheta = 25*eye(dim.n_theta); %25 (broad) variance identity matrix -- leads to sd of 5, and 1/(1+exp(10)) approaches zero

%phi is K (1), lambda (2), nu (3), rho (4), epsilon (5)
priors.muPhi = [0, ... %K (intercept) is ~N(0,1) transformed into ~unif(0,4000). prior mu of 0 corresponds to middle of distribution = 2000
    0, ... %lambda (autocorr) is ~N(0, 25) with exponential transform. So prior mu of 0 means lambda = 0.5
    0, ... %nu (go for gold) is ~N(0, 25) with exponential transform and rescaled to max = 100. Prior mu of 0 means nu = 50
    -5, ... %rho (mean fast/slow) is ~N(0, 1), transformed into ~gamma(2,2)*400, so a value of -5 corresponds to a rho of 0.5 (close to zero since there is no reason to expect mean diff a priori)
    -3, ... %epsilon (sd fast-slow explore) is ~N(0, 1), transformed into ~gamma(2,2)*1600, so a value of -3 is an epsilon of 169 (weak explore)
];

priors.SigmaPhi = diag([1, 25, 25, 1, 1]); %initial variance of the phi parameters

priors.muX0 = [mean(rts), ... %initial value of best RT 
    1.01, ... %a_fast (initial beta distribution hyperparameters)
    1.01, ... %b_fast
    1.01, ... %a_slow
    1.01, ... %b_slow
    0, ... %V (value)
    0, ... %Go
    0, ... %NoGo
    mean(rts), ... %initial RTlocavg
    0]; %sign of exploration influence

priors.SigmaX0 = zeros(dim.n); %0 variance in initial states

options.priors = priors;

%% set up models within evolution/observation Fx
options.inF.maxAlpha = 5; %maximum alpha for Go and NoGo updates
options.inG.maxNu = 100; %maximum nu for go for gold
options.inG.rhoMultiply = 400; %leads to a max rho around 10000
options.inG.epsilonMultiply = 1600; %leads to a max epsilon around 58000
options.inF.priors = priors; %copy priors into inF for parameter transformation
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)
options.inG.maxRT = 4000; %maximum possible RT (used for uniform distribution)
options.inG.meanRT = mean(rts);
%options.inG.multinomial = multinomial;


options.TolFun = 1e-8;
options.GnTolFun = 1e-8;
options.verbose=1;
options.DisplayWin=1;


rtsbyrun = reshape(rts, n_t/n_runs, n_runs)'; %convert to 8 x 50
runsum = [0, cumsum(options.multisession.split)];
rtslag = NaN(n_runs, n_t/n_runs);
for i = 1:(length(runsum)-1)
    rtslag(i, :) = rts([runsum(i)+1, (runsum(i)+1):(runsum(i+1)-1)]);
end
rtslag = reshape(rtslag', n_t, 1); %flatten into vector row-wise

rewardsbyrun = reshape(rewards, n_t/n_runs, n_runs)'; %convert to 8 x 50
rew_max = NaN(n_runs, n_t/n_runs);
for i = 1:n_runs
    best_rew=0;
    for t = 1:size(rewardsbyrun, 2)
        if rewardsbyrun(i, t) > best_rew, best_rew = rewardsbyrun(i, t); end
        rew_max(i, t) = best_rew;
    end
end
rew_max = reshape(rew_max', n_t, 1); %flatten into vector row-wise

rew_std = NaN(n_runs, n_t/n_runs);
for i = 1:n_runs
    for t = 1:size(rewardsbyrun, 2)
        rew_std(i, t) = std(rewardsbyrun(i, 1:t));
    end
end
rew_std = reshape(rew_std', n_t, 1); %flatten into vector row-wise

u  = [rts'; rtslag'; rewards'; rew_max'; rew_std'];

[posterior,out] = VBA_NLStateSpaceModel(y,u,@clock_tc_evolution,@clock_tc_observation,dim,options);

%% save output figure
h = figure(1);
%savefig(h,sprintf('%d_%s_multinomial%d_multisession%d_fixedParams%d',id,model,multinomial,multisession,fixed_params_across_runs))
savefig(h,sprintf('%d_tc_vba',id))
save(sprintf('%d_tc_vba_fit', id), 'posterior', 'out');
