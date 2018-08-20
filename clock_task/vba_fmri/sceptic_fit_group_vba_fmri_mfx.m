%loads in subjects' data and fits SCEPTIC models using VBA;

close all;
clear;
%curpath = fileparts(mfilename('fullpath'));

%addpath(genpath('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task'));
% addpath('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task'); %has glob.m and setup_rbf.m
% addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox'));

ncpus=getenv('matlab_cpus');

os = computer;
[~, me] = system('whoami');
me = strtrim(me);
is_alex=strcmp(me,'Alex')==1;


if is_alex
  addpath('~/code/temporal_instrumental_agent/clock_task'); %has glob.m and setup_rbf.m
  addpath(genpath('~/code/VBA-toolbox'));
  cd('~/code/temporal_instrumental_agent/clock_task/subjects')
  behavfiles = glob('*.csv');
  group_dir = '~/code/temporal_instrumental_agent/clock_task/subjects';
  ncpus=4;
  fprintf('defaulting to 4 cpus on old iMac \n');
else
  behavfiles = glob('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/subjects/*.csv');
  % behavfiles = glob('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/subjects/SPECC/*.csv');
  % behavfiles = glob('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/subjects/*.csv');

  group_dir = '/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/subjects';
  addpath('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task'); %has glob.m and setup_rbf.m
  addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox'));
  if strcmpi(ncpus, '')
    ncpus=40;
    fprintf('defaulting to 40 cpus because matlab_cpus not set\n');        
  else
    ncpus=str2double(ncpus);
  end    
end

%% set parameters
nbasis = 24;
multinomial = 1;
multisession = 0;
fixed_params_across_runs = 1;
fit_propspread = 0;
n_steps = 40; %number of time steps (size of multinomial vector)
u_aversion=0; %not used in fmri code at present
graphics = 0; %don't display interactive fitting

% get ID list
ids = cell(length(behavfiles),1);

%% main loop
<<<<<<< HEAD
if is_alex
  poolobj = gcp('nocreate');
  if isempty(poolobj)
    poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
  end
else
  ncpus=getenv('matlab_cpus');
  if strcmpi(ncpus, '')
    ncpus=40;
    fprintf('defaulting to 40 cpus because matlab_cpus not set\n');
  else
    ncpus=str2double(ncpus);
  end

  poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
  
end

model = 'fixed'; % will run to get value and prediction errors.

if strcmp(model,'fixed')
  f_fname = @h_sceptic_fixed_fmri;
  hidden_states=2;
elseif strcmp(model,'decay')
  f_fname = @h_sceptic_fixed_decay_fmri;
  hidden_states=3;
end

y_all = cell(length(behavfiles), 1);
u_all = cell(length(behavfiles), 1);
options_all = cell(length(behavfiles), 1);

for sub = 1:length(behavfiles)
    [~, str] = fileparts(behavfiles{sub});
    ids{sub} = regexp(str,'(?<=fMRIEmoClock_)[\d_]+(?=_tc)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

    fprintf('Loading subject %d id: %s \r',sub, char(ids{sub}));
    [data, y, u] = sceptic_get_data(behavfiles{sub}, n_steps);
    [options, dim] = sceptic_get_options(data, nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps, u_aversion, graphics);
    
    % populate data structures for VBA_MFX
    y_all{sub} = y;
    u_all{sub} = u;
    options_all{sub} = options;
end

%options for MFX
options_group.TolFun=2e-2;
options_group.MaxIter=50;
options_group.DisplayWin=0;
options_group.verbose=1;

priors_group.muPhi = 0; %temperature -- exp(phi(1))
priors_group.SigmaPhi = 10; %variance on temperature (before exponential transform)
priors_group.muTheta = [0; 0]; %learning rate (alpha), selective maintenance (gamma) -- before logistic transform
priors_group.SigmaTheta =  1e1*eye(10); %variance of 10 on alpha and gamma
priors_group.muX0 = zeros(nbasis*hidden_states,1); %have PE and decay as tag-along states
priors_group.SigmaX0 = zeros(nbasis*hidden_states, nbasis*hidden_states); %have PE and decay as tag-along states
priors_group.a_vX0 = repmat(Inf, [1, nbasis*hidden_states]); %use infinite precision prior on gamma for X0 to treat as fixed (a = Inf; b = 0)
priors_group.b_vX0 = repmat(0, [1, nbasis*hidden_states]);

[p_sub, o_sub, p_group, o_group] = VBA_MFX_parallel(y_all, u_all, f_fname, @g_sceptic, dim, options_all, priors_group, options_group);

%[p_sub, o_sub, p_group, o_group] = VBA_MFX(y_all, u_all, @h_sceptic_fixed_decay_fmri, @g_sceptic, dim, options_all, priors_group, options_group);

if is_alex
    cd ~/'Box Sync'/skinner/projects_analyses/SCEPTIC/mfx_analyses/
    if ~=exist(model,'dir')
        mkdir(model)
        cd(model)
    end
end

%too huge to save into one .mat file
save('vba_mfx_results_psub.mat', 'p_sub', '-v7.3');
save('vba_mfx_results_pgroup.mat', 'p_group', '-v7.3');
save('vba_mfx_results_ogroup.mat', 'o_group', '-v7.3');
save('vba_mfx_results_osub.mat', 'o_sub', '-v7.3');
save('vba_mfx_results_settings.mat', 'priors_group', 'options_group', 'y_all', 'u_all', 'ids', '-v7.3');

delete(poolobj);
