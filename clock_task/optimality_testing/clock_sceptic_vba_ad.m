close all
clear variables
n_basis = 8; %24%
range_RT = 400;
n_steps = 40;
n_t = 400;
n_runs = 8;
trialsToFit = 1:n_t;

%% fit as multiple runs
multisession = 1;
% fix parameters across runs
fixed = 1;

%u is 2 x ntrials where first row is rt and second row is reward
data = readtable('/Users/localadmin/code/clock_smoothoperator/clock_task/subjects/fMRIEmoClock_10637_tc_tcExport.csv','Delimiter',',','ReadVariableNames',true);
% data = readtable('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/fMRIEmoClock_10637_tc_tcExport.csv','Delimiter',',','ReadVariableNames',true);

dim = struct('n',n_basis,'n_theta',1,'n_phi',1,'p',n_steps, 'n_t', n_t);

if multisession
    options.multisession.split = repmat(n_t/n_runs,1,n_runs); % two sessions of 120 datapoints each
    %% fix parameters
    if fixed
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
%         
        % allow unique initial values for each trustee?
%         options.multisession.fixed.X0 = 'all';
    end
    
    
end

priors.muTheta = zeros(dim.n_theta,1);
priors.muPhi = zeros(dim.n_phi,1); % exp tranform
priors.muX0 = zeros(dim.n,1);
priors.SigmaPhi = 1e1*eye(dim.n_phi); %identifiability
priors.SigmaTheta = 1e1*eye(dim.n_theta);
priors.SigmaX0 = zeros(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
options.inF.nbasis = n_basis;
options.inF.ntimesteps = n_steps;
options.inG.ntimesteps = n_steps;
options.TolFun = 1e-8;
options.GnTolFun = 1e-8;

[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, .08);


rtrnd = round(data{trialsToFit,'rt'}*0.1*n_steps/range_RT)';
rtbin = zeros(n_steps, length(trialsToFit));
for i = trialsToFit
   rtbin(rtrnd(i), i) = 1;
end

options.sources(1) = struct('out',1:n_steps,'type',2);

u = [(data{trialsToFit, 'rt'}*0.1*n_steps/range_RT)'; data{trialsToFit, 'score'}'];

[posterior,out] = VBA_NLStateSpaceModel(rtbin,u,@clock_sceptic_evolution,@clock_sceptic_observation,dim,options);

