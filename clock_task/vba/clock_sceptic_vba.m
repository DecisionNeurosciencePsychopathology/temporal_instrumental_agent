function [posterior,out] = clock_sceptic_vba(id,model,n_basis,multinomial,multisession,fixed_params_across_runs)

%% fits SCEPTIC model to Clock Task subject data using VBA toolbox

% example call:
% [posterior,out]=clock_sceptic_vba(10638,'fixed',4,1,1,1)

% id:           5-digit subject id in Michael Hallquist's BPD study
% only works with 'fixed' (fixed learning rate SCEPTIC) so far
% n_basis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)

%%
close all

global rew_rng_state
rew_rng_seed = 99;

% n_basis = 8; %24%
options.inF.fit_nbasis = 0;
range_RT = 400;
n_steps = 40;
n_t = 400;
n_runs = 8;
trialsToFit = 1:n_t;

%% fit as multiple runs
% multisession = 1;
% fix parameters across runs
% fixed_params_across_runs = 1;

%% u is 2 x ntrials where first row is rt and second row is reward
% If we can't find the path have the user select it.
try
data = readtable(sprintf('c:/kod/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_%d_tc_tcExport.csv', id),'Delimiter',',','ReadVariableNames',true);
% data = readtable('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_10637_tc_tcExport.csv','Delimiter',',','ReadVariableNames',true);
catch
    disp('Can''t read subject data!\n\n')
    [fname,path_name]=uigetfile('*.csv', 'Select a subjects data file');
    data = readtable([path_name fname],'Delimiter',',','ReadVariableNames',true);
end


%% set up models within evolution/observation Fx
%Note: we might need to add option.inF.model to make the kalman models
%easier to deal with...
options.inF.nbasis = n_basis;
options.inF.ntimesteps = n_steps;
options.inG.ntimesteps = n_steps;
options.inG.multinomial = multinomial;
options.inG.nbasis = n_basis;

%%
options.TolFun = 1e-8;
options.GnTolFun = 1e-8;
options.verbose=1;
options.DisplayWin=1;

%% set up basis
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, .08);


%Set up sigma noise for every point in u or hidden state?
    rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
    rew_rng_state=rng;
    [~,idx] = unique(data.run);
    conditions=data.rewFunc(idx);
    sigma_noise = [];
    bin_size = n_t/n_runs;
    for i = 1:length(conditions)
        sigma_noise = [sigma_noise repmat(std(arrayfun(@(x) RewFunction(x*10, conditions(i), 0), options.inF.tvec))^2, 1, bin_size)];
    end
    sigma_noise = mean(sigma_noise);
options.inF.sigma_noise = sigma_noise;


%% split into conditions/runs
if multisession %improves fits moderately
    options.multisession.split = repmat(n_t/n_runs,1,n_runs); % two sessions of 120 datapoints each
    %% fix parameters
    if fixed_params_across_runs
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        %
        % allow unique initial values for each run?
        %                  options.multisession.fixed.X0 = 'all';
    end
    
end



%Determine which evolution funciton to use
switch model 
    case 'fixed'
        h_name = @h_sceptic_fixed;
        hidden_variables = 1; %tracks only value
        priors.muX0 = zeros(hidden_variables*n_basis,1);
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        
    case 'kalman_softmax'
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        h_name = @h_sceptic_kalman_softmax;
    case 'kalman_processnoise'
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        h_name = @h_sceptic_kalman_processnoise;
    case 'kalman_uv_sum'
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis);
        h_name = @h_sceptic_kalman_uv_sum;
end
% if strcmpi(model,'fixed')
%     h_name = @h_sceptic_fixed;
% elseif strcmpi(model,'kalman_softmax')
%     h_name = @h_sceptic_kalman_softmax; % Jon will continue here
% end
g_name = @g_sceptic;



if multinomial
    rtrnd = round(data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
    dim = struct('n',hidden_variables*n_basis,'n_theta',1,'n_phi',1,'p',n_steps);
    options.sources(1) = struct('out',1:n_steps,'type',2);
    
    %% compute multinomial response -- renamed 'y' here instead of 'rtbin'
    y = zeros(n_steps, length(trialsToFit));
    for i = trialsToFit
        y(rtrnd(i), i) = 1;
    end
    priors.a_alpha = Inf;   % infinite precision prior
    priors.b_alpha = 0;
    priors.a_sigma = 1;     % Jeffrey's prior
    priors.b_sigma = 1;     % Jeffrey's prior
    options.binomial = 1;
    priors.SigmaPhi = 1e1*eye(dim.n_phi);
else
    dim = struct('n',hidden_variables*n_basis,'n_theta',1,'n_phi',1, 'n_t', n_t);
    y = (data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
    priors.a_sigma = 1;     % Jeffrey's prior
    priors.b_sigma = 1;     % Jeffrey's prior
    priors.SigmaPhi = 0*eye(dim.n_phi);
    options.binomial = 0;
    options.sources(1) = struct('out',1,'type',0);
    
    
end
%
% if options.inF.fit_nbasis
%     dim = struct('n',n_basis,'n_theta',2,'n_phi',1,'p',n_steps);
% priors.muTheta = [0 8];
% priors.muPhi = zeros(dim.n_phi,1); % exp tranform
% priors.muX0 = zeros(dim.n,1);
% priors.SigmaPhi = 1e1*eye(dim.n_phi);
% priors.SigmaTheta = 1e1*eye(dim.n_theta);
% options.inF.priordist_theta2 = makedist('Normal',priors.muTheta(2), unique(max(priors.SigmaTheta)));
% options.inF.maxbasis = 24;
% options.inF.muTheta2 = priors.muTheta(2);
% options.inF.SigmaTheta2 = unique(max(priors.SigmaTheta));
% else
% priors.muTheta = zeros(dim.n_theta,1);
% priors.muPhi = zeros(dim.n_phi,1); % exp tranform
% priors.muX0 = zeros(dim.n,1);
% priors.SigmaPhi = 1e1*eye(dim.n_phi);
% priors.SigmaTheta = 1e1*eye(dim.n_theta);

% end
%% skip first trial
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;

%% priors
priors.muTheta = zeros(dim.n_theta,1);
priors.muPhi = zeros(dim.n_phi,1); % exp tranform
priors.SigmaTheta = 1e1*eye(dim.n_theta);
options.priors = priors;


%u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'];
u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'];

[posterior,out] = VBA_NLStateSpaceModel(y,u,h_name,g_name,dim,options);

%% save output figure
h = figure(1);
savefig(h,sprintf('%d_%s_multinomial%d_multisession%d_fixedParams%d',id,model,multinomial,multisession,fixed_params_across_runs))