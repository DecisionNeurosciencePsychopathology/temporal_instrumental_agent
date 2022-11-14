function [options, dim] = sceptic_get_vba_options(data, so)
% this function sets up the options and dim structures for VBA
% so should be setup upstream by sceptic_validate_options

options=[];
priors=[];

if ~so.graphics
  options.DisplayWin = 0; %whether to display graphics during fitting
  options.GnFigs = 0;
end

% u is 2 x ntrials where first row is rt and second row is reward

% copy so to inF and inG as starting point
options.inF = so;
options.inG = so;

options.inF.fit_nbasis = 0; %do not try different numbers of basis functions in fitting

% convergence settings
options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=0; %don't show single subject fitting process

%uses max prop spread parameter to obtain refspread in case where fit_propspread = 0;
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread,options.inG.gaussmat_highres, options.inF.gaussmat_trunc_highres] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, options.inF.max_prop_spread);

options.inF.gaussmat = options.inG.gaussmat; %copy (non-truncated) gaussmat into inF

n_t = size(data, 1); %number of rows
n_runs = n_t/so.trials_per_run; %determine number of runs
options.inF.n_runs = n_runs; %used downstream

%% split into conditions/runs
if so.multisession
  options.multisession.split = repmat(n_t/n_runs,1,n_runs);
  
  % fix parameters
  if fixed_params_across_runs
    options.multisession.fixed.theta = 'all';
    options.multisession.fixed.phi = 'all';
    
    % allow unique initial values for each run?
    options.multisession.fixed.X0 = 'all';
  end
end

%% skip evolution on first trial (copy initial states)
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;

%% specify dimensions of data to be fit
dim = struct('n', so.hidden_states * so.nbasis, ...
  'n_theta', so.n_theta + so.fit_propspread, ...
  'n_phi', so.n_phi, ...
  'p', so.ntimesteps, ...
  'n_t', n_t);

%%populate priors
priors = sceptic_get_priors(dim, so);

options.priors = priors;
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)

if so.multinomial
  options.sources(1) = struct('out', 1:so.ntimesteps, 'type', 2); %multinomial fitting
end

end
