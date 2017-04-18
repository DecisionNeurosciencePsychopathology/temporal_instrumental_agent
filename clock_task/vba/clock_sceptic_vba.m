function [posterior,out] = clock_sceptic_vba(s,id,model, data_file)

%% fits SCEPTIC model to Clock Task subject data using VBA toolbox
% example call:
% [posterior,out] = clock_sceptic_vba(s,10637,'fixed, 'subjects/fMRIEmoClock_10637_tc_tcExport.csv')
% s:            Configuration struct (see create_scecptic_configuration_struct)
% id:           5-digit subject id in Michael Hallquist's BPD study
% model:        Current model to be fit
% data_file:    String of subject's path to processed data
% n_basis:      Number of basis functions
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
% fixed_params_across_runs: self-explanatory
% fit_propspread: makes temporal generalization within the eligibility trace a free parameter
% n_steps:      number of time bins
% u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum
% saveresults:  self-explanatory
% graphics:     Plot VBA graphics or  not
% results_dir:  Where to store saved results
% cstruct:      Condition struct
% task_name:    Which version of clock task is being used
% range_RT:     Range of maximum RT subject could have

%% Set up
%Close all existing figures
close all

%Must have the config structure, id, and model to run.
try
    s; %Config struct
    id; %Subj id
    model; %Current model
    data_file; %Subj data file path
catch
    error('Main function arguments not found! Exiting...')
end

%Set all parameters into local variables from input struct
try n_basis=s.n_basis; catch, n_basis = 24; end
try multinomial=s.multinomial; catch, multinomial = 1; end
try multisession=s.multisession; catch, s.multisession = 0; end
try fixed_params_across_runs=s.fixed_params_across_runs; catch, fixed_params_across_runs = 0; end
try fit_propspread=s.fit_propspread; catch, fit_propspread = 0; end
try n_steps=s.n_steps; catch, n_steps = 50; end
try u_aversion=s.u_aversion; catch, u_aversion = 1; end
try saveresults=s.saveresults; catch, saveresults = 0; end
try graphics=s.graphics; catch, graphics = 0; end
try results_dir=s.results_dir; catch, results_dir = pwd; end
try cstruct=s.cstruct; catch, cstruct = []; end
try task_name=s.task_name; catch, task_name = 'hallquist_clock'; end
try range_RT=s.range_RT; catch, range_RT = 400; end


%Set reward seeds
global rew_rng_state no_gamma
rew_rng_seed = 99;

usecstruct=0;
if isstruct(cstruct)
    cond = cstruct.name;
    usecstruct=1;
end

%Graphics handle
if ~graphics
    options.DisplayWin = 0;
    options.GnFigs = 0;
end

%% set up dim defaults
n_theta = 1;
n_phi = 1;

%% no choice autocorrelation by default
%Autocorrelation choices
%1) 'none'                      No autocorrelation
%2) 'softmax_multitrial'        Implements choice autocorrelation as in Schoenberg et al. 2007 without temporal generalization
%3) 'softmax_multitrial_smooth' Implements choice autocorrelation as in Schoenberg et al. 2007 with temporal generalization controlled by an additional temporal smoothing parameter iota
%4) 'exponential'
%5) 'choice_tbf'
options.inG.autocorrelation = 'none';
options.inF.autocorrelation = options.inG.autocorrelation;

%Entropy
options.inF.entropy = 0; %If we want to track entropy per trial
track_entropy=options.inF.entropy;
options.inF.H_threshold = 0.01;

%If we want to use the elig update variant in choice rule
options.inF.total_pe=0;

%If we want to track delta for regressors
options.inF.track_pe = 1;

%% FILE I/O

% Import the subject's data
data = readtable(data_file,'Delimiter',',','ReadVariableNames',true);
n_t = size(data,1); %Number of trials
n_runs = n_t/50; %Number of runs
trialsToFit = 1:n_t;

%% set up models within evolution/observation Fx
options.inF.fit_nbasis = 0;
options.inF.fit_propspread = fit_propspread;
options.inF.nbasis = n_basis;
options.inF.ntimesteps = n_steps;
options.inG.ntimesteps = n_steps;
options.inG.multinomial = multinomial;
options.inG.nbasis = n_basis;
options.inG.maxRT = range_RT;

options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=1;

%% set up kalman defaults
options.inF.kalman.kalman_softmax = 0;
options.inF.kalman.kalman_uv_sum = 0;
options.inF.kalman.fixed_uv = 0;

%% set up basis
fixed_prop_spread = .0125;
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, fixed_prop_spread);
options.inG.sig_spread = options.inF.sig_spread;

rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
rew_rng_state=rng;
sigma_noise = [];
if usecstruct
    %taking std over all timesteps and possible draws here. This is in contrast to approach below where you get one
    %sample from the contingency as drawn by RewFunction. This is much faster than arrayfun below.
    sigma_noise = repmat(std(cstruct.lookup(:))^2, 1, n_basis);
else
    [~,idx] = unique(data.run);
    conditions=data.rewFunc(idx);
    
    run_length = n_t/n_runs;
    for i = 1:length(conditions)
        sigma_noise = [sigma_noise repmat(std(arrayfun(@(x) RewFunction(x*100, conditions(i), 0), options.inF.tvec))^2, 1, run_length)];
    end
end

sigma_noise = mean(sigma_noise);
options.inF.sigma_noise = sigma_noise;
options.inF.gaussmat = options.inG.gaussmat;

%% split into conditions/runs
if multisession %improves fits moderately
    options.multisession.split = repmat(n_t/n_runs,1,n_runs); % two sessions of 120 datapoints each
    %fix parameters
    if fixed_params_across_runs
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        options.multisession.fixed.X0 = 'all';
    end
end

%% Determine which evolution funciton to use
options.inF.kalman.(model)=1; %Declare which model to use if kalman

switch model
    case 'fixed'
        h_name = @h_sceptic_fixed;
        hidden_variables = 1; %tracks only value
        priors.muX0 = zeros(hidden_variables*n_basis,1);
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
    
    case 'fixed_uv'
        h_name = @h_sceptic_kalman;
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        n_theta = 2; %tau alpha
        options.inF.u_aversion = u_aversion;
        options.inG.u_aversion = u_aversion;
        
    case 'fixed_decay'
        h_name = @h_sceptic_fixed_decay;
        hidden_variables = 1; %tracks only value
        priors.muX0 = zeros(hidden_variables*n_basis,1);
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        n_theta = 2; %learning rate and decay outside of the eligibility trace
        
    case 'kalman_softmax'
        h_name = @h_sceptic_kalman;
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        %Prop_spread is the only variable in this model
        if fit_propspread
            n_theta = 0;
        end

    case 'kalman_uv_sum'
        h_name = @h_sceptic_kalman;
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        options.inF.u_aversion = u_aversion;
        options.inG.u_aversion = u_aversion;
        
    otherwise
        disp('The model you have entered does not match any of the default names, check spelling!');
        return
        
end

%Adjust initial arrays if tracking pe
if options.inF.track_pe == 1
    track_pe = 1;
    hidden_variables = hidden_variables + 1;
    priors.muX0 = [priors.muX0; zeros(n_basis,1)];
    priors.SigmaX0 = zeros(hidden_variables*n_basis);
else
    track_pe = 0;
end

%Add in the lambda parameter if using autocorrelation
if strcmp(options.inG.autocorrelation,'exponential') || strcmp(options.inG.autocorrelation,'softmax_multitrial')
    n_phi = n_phi + 2;
elseif strcmp(options.inG.autocorrelation,'softmax_multitrial_smooth')
    n_phi = n_phi + 3;
elseif strcmp(options.inG.autocorrelation,'choice_tbf') %Modifiy the models if we want to add the choice basis autocorrelation
    n_phi = n_phi + 1;
    n_theta = n_theta+1; %Add choice decay parameter
    hidden_variables = hidden_variables + 1; %Add in choice basis
    priors.muX0 = [priors.muX0; zeros(n_basis,1)];
    priors.SigmaX0 = zeros(hidden_variables*n_basis);
elseif options.inF.entropy == 1
    track_entropy=2; %Add 2 hidden states to track entropy and max value pre-update
    priors.muX0 = [priors.muX0; zeros(track_entropy,1)];
    priors.SigmaX0 = zeros(hidden_variables*n_basis+track_entropy);
end

%Finalize hidden variables
options.inF.hidden_state = hidden_variables;

%Map the necessary options from F to G
options.inG.hidden_state = options.inF.hidden_state;
options.inG.kalman = options.inF.kalman;

if multinomial
    rtrnd = round(data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
    rtrnd(rtrnd==0)=1;
    dim = struct('n',hidden_variables*n_basis+track_entropy,'n_theta',n_theta+fit_propspread,'n_phi',n_phi,'p',n_steps);
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
    u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'; trialsToFit];
    u = [zeros(size(u,1),1) u(:,1:end-1)];
    options.inG.rts = round(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)';
    % Observation function
    g_name = @g_sceptic;

else
    n_phi = 2; % [autocorrelation lambda and response bias/meanRT K] instead of temperature
    dim = struct('n',hidden_variables*n_basis,'n_theta',n_theta+fit_propspread,'n_phi',n_phi, 'n_t', n_t);
    y = (data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
    priors.a_sigma = 1;     % Jeffrey's prior
    priors.b_sigma = 1;     % Jeffrey's prior
    priors.muPhi = [0, 0];  % K, lambda
    %priors.SigmaPhi = diag([0,1]); % get rid of the K
    priors.SigmaPhi = diag([1,1]);
    options.binomial = 0;
    options.sources(1) = struct('out',1,'type',0);
    prev_rt = [0 y(1:end-1)];
    % Inputs
    u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'; prev_rt];
    u = [zeros(size(u,1),1) u(:,1:end-1)];
    
    % Observation function
    g_name = @g_sceptic_continuous;
    
end

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

%% Run VBA State Space Model
[posterior,out] = VBA_NLStateSpaceModel(y,u,h_name,g_name,dim,options);

%% save output
if saveresults
    save([results_dir sprintf(s.results_str, s.results_str_values, id, model)], 'posterior', 'out', '-v7.3');
end
