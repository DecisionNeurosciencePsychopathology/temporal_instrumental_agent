function [posterior,out] = clock_tc_klambda(id,showfig)
%% fits K + lambda model to Clock Task subject data using VBA toolbox
% example call:
% [posterior,out]=clock_sceptic_vba(10638,'fixed',4,1,1,1)
% id:           5-digit subject id in Michael Hallquist's BPD study
% only works with 'fixed' (fixed learning rate SCEPTIC) so far
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)

close all

n_t = 400; %400 trials
n_runs = 8;
trialsToFit = 1:n_t;

if nargin < 2, showfig = 1; end

if showfig==1
    options.DisplayWin=1;
else
    options.DisplayWin=0;
end

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
        % options.multisession.fixed.X0 = 'all';
    end
end

n_states = 0;

dim = struct('n',n_states,'n_theta',0,'n_phi',2, 'n_t', n_t);
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

%% priors
rtrescale = 1000; %put the RT-related hidden states on a similar scale as the others (i.e., divide by 1000)

% 20Jan2016: After discussing the issue of priors with Alex, we agreed to try setting the prior Gaussian means closer to 0 to produce an 
% uninformative prior. This leads to the priors being placed in the middle of the parameter space in each case, regardless of parameter
% transformation. So, for example, for the epsilon parameter, which varies from ~0..~

%phi is K (1), lambda (2)
% For uniform and gamma distribution transformations using the iCDF method, use the standard normal distribution to allow for precompiled
% fastnormcdf function to work.

priors.muPhi = [0, ... %K (intercept) is ~N(0, 1) transformed into ~unif(0,4000). prior mu of 0 corresponds to middle of distribution = 2000
    0 ]; %lambda (autocorr) is ~N(0, 15) with exponential transform. So prior mu of 0 means lambda = 0.5

priors.SigmaPhi = diag([1, 15]); %initial variance of the phi parameters

options.priors = priors;

%% set up models within evolution/observation Fx
options.inF.priors = priors; %copy priors into inF for parameter transformation
options.inG.maxNu = 10; %maximum nu for go for gold (tends to be between 0 and 1, though)
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)
options.inG.maxRT = 4000; %maximum possible RT (used for uniform distribution)
options.inG.meanRT = mean(rts);
options.inG.RTrescale = rtrescale; %scale bestRT hidden state into seconds (for similarity across hidden states)

options.TolFun = 1e-8;
options.GnTolFun = 1e-8;
options.verbose=1;

rtsbyrun = reshape(rts, n_t/n_runs, n_runs)'; %convert to 8 x 50
runsum = [0, cumsum(options.multisession.split)];
rtslag = NaN(n_runs, n_t/n_runs);
for i = 1:(length(runsum)-1)
    rtslag(i, :) = rts([runsum(i)+1, (runsum(i)+1):(runsum(i+1)-1)]);
end
rtslag = reshape(rtslag', n_t, 1); %flatten into vector row-wise

u  = [rts'; rtslag'];

%Gaussian -> Uniform transformation using the inverse CDF method is very slow
%Is it better to precompute the Gaussian CDF using the priors, then use a cubic spline interpolation lookup?
%This yields very small error (on the order of 1e-17).
%it turns out that the repeated calls to interpn are slower than repeated cdf evaluations
%options.inG.stdnormcdf = cdf('Normal', -10:.0001:10, 0, 1);
%options.inG.stdnormgrid = -10:.0001:10; %sampling grid

tic;
[posterior,out] = VBA_NLStateSpaceModel(y,u,@clock_klambda_evolution,@clock_klambda_observation,dim,options);
elapsed=toc;
fprintf('Model converged in %.2f seconds\n', elapsed);

posterior.transformed.K = unifinv(fastnormcdf(posterior.muPhi(1)), 0, options.inG.maxRT);
posterior.transformed.lambda = 1 / (1+exp(-posterior.muPhi(2)));

%% save output figure
if (showfig==1)
    h = figure(1);
    savefig(h,sprintf('results/%d_klambda_vba',id));
    %savefig(h,sprintf('%d_%s_multinomial%d_multisession%d_fixedParams%d',id,model,multinomial,multisession,fixed_params_across_runs))
end
save(sprintf('results/%d_klambda_vba_fit', id), 'posterior', 'out');
