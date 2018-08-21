function [options, dim] = sceptic_get_vba_options(data, so)
  %% this function sets up the options and dim structures for VBA
  
  so=sceptic_validate_options(so); %should be handled upstream, but just in case
  
options=[];
priors=[];

%% set up dim defaults
n_phi = 1; %temperature

if ~so.graphics
  options.DisplayWin = 0;
  options.GnFigs = 0;
end

%% u is 2 x ntrials where first row is rt and second row is reward
% If we can't find the path have the user select it.

options.inF.fit_nbasis = 0; %do not try different numbers of basis functions in fitting
options.inF.fit_propspread = so.fit_propspread;

%% initialize key setting within evolution/observation Fx inputs
options.inF.nbasis = so.nbasis;
options.inF.ntimesteps = so.ntimesteps;
options.inG.ntimesteps = so.ntimesteps;
options.inG.multinomial = so.multinomial;
options.inG.nbasis = so.nbasis;

%options.inG.maxRT = range_RT; %not used at present

%%convergence settings
options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=0; %don't show single subject fitting process

options.DisplayWin = so.graphics; %whether to display graphics during fitting

options.inF.max_prop_spread = so.max_prop_spread;

%uses max prop spread parameter to obtain refspread in case where fit_propspread = 0;
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, options.inF.max_prop_spread);

options.inF.gaussmat = options.inG.gaussmat; %copy (non-truncated) gaussmat into inF

n_t = size(data,1); %number of rows
n_runs = n_t/so.trials_per_run; %determine number of runs
options.inF.n_runs = n_runs; %used downstream

%% split into conditions/runs
if so.multisession
  options.multisession.split = repmat(n_t/n_runs,1,n_runs);

  %% fix parameters
  if fixed_params_across_runs
    options.multisession.fixed.theta = 'all';
    options.multisession.fixed.phi = 'all';
    
    %% allow unique initial values for each run?
    options.multisession.fixed.X0 = 'all';
  end
end

% setup number of hidden variables
options.inF.hidden_states = so.hidden_states; %used to index hidden state vector inside evolution and observation functions

%Map the necessary options from F to G
options.inG.hidden_states = options.inF.hidden_states;

%% skip first trial
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;

%% specify dimensions of data to be fit
dim = struct('n', so.hidden_states * so.nbasis, ...
	     'n_theta', so.n_theta + so.fit_propspread, ...
	     'n_phi', n_phi, ...
	     'p', so.ntimesteps, ...
	     'n_t', n_t);

%%populate priors
priors = sceptic_get_priors(dim);

options.priors = priors;
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)

if so.multinomial
  options.sources(1) = struct('out', 1:so.ntimesteps, 'type', 2);
  options.binomial = 1; %multinomial fitting
end

end