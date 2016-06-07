function [posterior,out] = clock_sceptic_vba_fmri(id,model,n_basis,multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,u_aversion, saveresults, graphics)
%% fits SCEPTIC model to single-subject Clock Task subject data using VBA toolbox
% example call:
% [posterior,out]=clock_sceptic_vba(10638,'modelname',nbasis,multinomial,multisession,fixed_params_across_runs,fit_propsrpead)
% id:           5-digit subject id in Michael Hallquist's BPD study
% only works with 'fixed' (fixed learning rate SCEPTIC) so far
% n_basis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
% fixed_params_across_runs -- self-explanatory
% fit_propspread -- makes temporal generalization within the eligibility trace a free parameter
% n_steps:      number of time bins
% u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum

close all

%% uncertainty aversion for UV_sum
if nargin < 7, fit_propspread = 0; end
if nargin < 9, u_aversion = 0; end
if nargin < 10, saveresults = 1; end
if nargin < 11, graphics = 0; end

global rew_rng_state 
rew_rng_seed = 99;

if ~graphics
    options.DisplayWin = 0;
    options.GnFigs = 0;
end
%% set up dim defaults
n_theta = 1;
n_phi = 1;

%% fit as multiple runs
% multisession = 1;
% fix parameters across runs
% fixed_params_across_runs = 1;

basedir = '/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/subjects';
results_dir = '/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/vba_fmri/vba_out';
data = readtable(sprintf('%s/fMRIEmoClock_%d_tc_tcExport.csv', basedir, id),'Delimiter',',','ReadVariableNames',true);

%% u is 2 x ntrials where first row is rt and second row is reward
% If we can't find the path have the user select it.

options.inF.fit_nbasis = 0; %do not try different numbers of basis functions in fitting
range_RT = 400;
% n_steps = 4000;
n_t = size(data,1);
n_runs = n_t/50; %50 trials per run
trialsToFit = 1:n_t;
options.inF.fit_propspread = fit_propspread;

%% set up models within evolution/observation Fx
%Note: we might need to add option.inF.model to make the kalman models
%easier to deal with...
options.inF.nbasis = n_basis;
options.inF.ntimesteps = n_steps;
options.inG.ntimesteps = n_steps;
options.inG.multinomial = multinomial;
options.inG.nbasis = n_basis;
options.inG.maxRT = range_RT;
%%
options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=1;
% options.DisplayWin=1;

%[c, sig, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(40, 16, .08);
%save('sceptic_fmri_basis_setup.mat', 'c', 'sig', 'tvec', 'sig_spread', 'gaussmat', 'gaussmat_trunc', 'refspread');

%% set up basis (fixed prop_spread of .08)?
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, .08);

%Set up sigma noise for every point in u or hidden state?
rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
rew_rng_state=rng;
[~,idx] = unique(data.run);
conditions=data.rewFunc(idx);
sigma_noise = zeros(length(conditions), 1);
run_length = n_t/n_runs;
for i = 1:length(conditions)
    sni = repmat(std(arrayfun(@(x) RewFunction(x*100, conditions(i), 0), options.inF.tvec))^2, 1, run_length);
    sigma_noise(i) = mean(sni);
end
sigma_noise = mean(sigma_noise);
options.inF.sigma_noise = sigma_noise;
options.inF.gaussmat = options.inG.gaussmat;

%% split into conditions/runs
if multisession %improves fits moderately
    options.multisession.split = repmat(n_t/n_runs,1,n_runs);
    %% fix parameters
    if fixed_params_across_runs
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        % allow unique initial values for each run?x
        options.multisession.fixed.X0 = 'all';
    end
    
end

h_name = @h_sceptic_fixed_decay_fmri;
hidden_variables = 3; %value, PE, decay
priors.muX0 = zeros(hidden_variables*n_basis,1);
priors.SigmaX0 = zeros(hidden_variables*n_basis);
n_theta = 2; %learning rate and decay outside of the eligibility trace

options.inF.hidden_state = hidden_variables;

%Map the necessary options from F to G
options.inG.hidden_state = options.inF.hidden_state;

rtrnd = round(data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
rtrnd(rtrnd==0)=1;
dim = struct('n',hidden_variables*n_basis,'n_theta',n_theta+fit_propspread,'n_phi',n_phi,'p',n_steps);
options.sources(1) = struct('out',1:n_steps,'type',2);

%% compute multinomial response -- renamed 'y' here instead of 'rtbin'
y = zeros(n_steps, length(trialsToFit));
for i = 1:length(trialsToFit)
    y(rtrnd(i), i) = 1;
end
priors.a_alpha = Inf;   % infinite precision prior
priors.b_alpha = 0;
priors.a_sigma = 1;     % Jeffrey's prior
priors.b_sigma = 1;     % Jeffrey's prior
options.binomial = 1;
priors.muPhi = zeros(dim.n_phi,1); % exp tranform
priors.SigmaPhi = 1e1*eye(dim.n_phi);
% Inputs

u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'];
u = [zeros(size(u,1),1) u(:,1:end-1)];

% Observation function
g_name = @g_sceptic;

%% skip first trial
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;

%% priors
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e1*eye(dim.n_theta); % lower the learning rate variance -- it tends to be low in the posterior
options.priors = priors;
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)

onsets=NaN(1,size(y,2));
for mm=1:size(y,2)
    pos = find(y(:,mm)==1);
    if isempty(pos), pos=NaN; end
    onsets(mm) = pos;
end

[posterior,out] = VBA_NLStateSpaceModel(y,u,h_name,g_name,dim,options);

%scepticrefit(posterior, out);

if saveresults
    %% save output figure
    % h = figure(1);
    % savefig(h,sprintf('results/%d_%s_multinomial%d_multisession%d_fixedParams%d',id,model,multinomial,multisession,fixed_params_across_runs))
    save(sprintf([results_dir, '/SHIFTED_U_CORRECT_%d_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d_sceptic_vba_fit'], id, model, multinomial,multisession,fixed_params_across_runs, u_aversion), 'posterior', 'out');
end

