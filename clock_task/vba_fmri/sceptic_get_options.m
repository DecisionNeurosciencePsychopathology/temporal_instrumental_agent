function [options, dim] = sceptic_get_options(data, n_basis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps, u_aversion, graphics)

if nargin < 2, n_basis=24; end
  if nargin < 3, multinomial=1; end
  if nargin < 4, multisession=0; end
  if nargin < 5, fixed_params_across_runs=1; end
  if nargin < 6, fit_propspread=1; end
  if nargin < 7, n_steps=40; end
  if nargin < 8, u_aversion=0; end
  if nargin < 9, graphics=0; end

options=[];
priors=[];

%% set up dim defaults
n_phi = 1; %temperature
n_theta = 2; %learning rate and selective maintenance parameters (decay outside of the eligibility trace)

if ~graphics
    options.DisplayWin = 0;
    options.GnFigs = 0;
end

%% u is 2 x ntrials where first row is rt and second row is reward
% If we can't find the path have the user select it.

options.inF.fit_nbasis = 0; %do not try different numbers of basis functions in fitting
options.inF.fit_propspread = fit_propspread;

%% set up models within evolution/observation Fx
options.inF.nbasis = n_basis;
options.inF.ntimesteps = n_steps;
options.inG.ntimesteps = n_steps;
options.inG.multinomial = multinomial;
options.inG.nbasis = n_basis;
%options.inG.maxRT = range_RT; %not used at present

%%convergence settings
options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=0; %don't show single subject fitting process

% options.DisplayWin=1;

%%for now, fix prop_spread to .0125.
options.inF.max_prop_spread = 0.0125;

%uses max prop spread parameter to obtain refspread in case where fit_propspread = 0;
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, options.inF.max_prop_spread);

options.inF.gaussmat = options.inG.gaussmat; %copy (non-truncated) gaussmat into inF

n_t = size(data,1); %number of rows
n_runs = n_t/50; %50 trials per run

%% split into conditions/runs
if multisession
  options.multisession.split = repmat(n_t/n_runs,1,n_runs);
  % fix parameters
  if fixed_params_across_runs
    options.multisession.fixed.theta = 'all';
    options.multisession.fixed.phi = 'all';
    % allow unique initial values for each run?
    options.multisession.fixed.X0 = 'all';
  end
end

% setup number of hidden variables
hidden_variables = 3; %value, PE, decay

options.inF.hidden_state = hidden_variables; %used to index hidden state vector inside evolution and observation functions

%Map the necessary options from F to G
options.inG.hidden_state = options.inF.hidden_state;

%% skip first trial
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;

%% specify dimensions of data to be fit
dim = struct('n',hidden_variables*n_basis,'n_theta',n_theta+fit_propspread,'n_phi',n_phi,'p',n_steps,'n_t',400);

%%populate priors
priors = sceptic_get_priors(dim);

options.priors = priors;
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)

options.sources(1) = struct('out',1:n_steps,'type',2);
options.binomial = 1; %multinomial fitting

end
