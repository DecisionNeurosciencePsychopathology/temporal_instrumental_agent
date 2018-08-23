%loads in subjects' data and fits SCEPTIC models using VBA;

close all;
clear;

sceptic_dataset = 'mmclock_fmri';
so.model = 'fixed_UV';
%curpath = fileparts(mfilename('fullpath'));

%addpath(genpath('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task'));
% addpath('/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task'); %has glob.m and setup_rbf.m
% addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox'));

ncpus=getenv('matlab_cpus');

os = computer;
[~, me] = system('whoami');
me = strtrim(me);
is_alex=strcmp(me,'Alex')==1;

%note that this function looks for 'sceptic_dataset' and 'sceptic_model'
%as environment variables so that this script can be scaled easily for batch processing
so = sceptic_validate_options(); %initialize and validate sceptic fitting settings

%% set environment and define file locations
if is_alex
  sceptic_repo='/Users/localadmin/code/temporal_instrumental_agent/clock_task';
  addpath(genpath('/Users/localadmin/code/VBA-toolbox')); %setup VBA
else
  sceptic_repo='/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task';
  addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox')); %setup VBA
end

addpath(sceptic_repo); %has glob.m and setup_rbf.m

group_dir=[sceptic_repo, '/subjects']; %not used for anything at the moment...
if strcmpi(so.dataset,'mmclock_meg')
  behavfiles = glob([sceptic_repo, '/subjects/mmclock_meg/*.csv']);
elseif strcmpi(so.dataset,'mmclock_fmri')
  behavfiles = glob([sceptic_repo, '/subjects/mmclock_fmri/*.csv']);
elseif strcmpit(so.dataset,'specc')
  behavfiles = glob([sceptic_repo, '/clock_task/subjects/SPECC/*.csv']);
end

% get ID list
ids = cell(length(behavfiles),1);

%% setup parallel parameters
if is_alex
  ncpus=4;
  fprintf('defaulting to 4 cpus on old iMac \n');

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

y_all = cell(length(behavfiles), 1);
u_all = cell(length(behavfiles), 1);
options_all = cell(length(behavfiles), 1);

for sub = 1:length(behavfiles)
    [~, str] = fileparts(behavfiles{sub});
    ids{sub} = regexp(str,'(?<=fMRIEmoClock_)[\d_]+(?=_tc)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

    fprintf('Loading subject %d id: %s \r', sub, char(ids{sub}));
    [data, y, u] = sceptic_get_data(behavfiles{sub}, so);
    [options, dim] = sceptic_get_vba_options(data, so);
    
    % populate data structures for VBA_MFX
    y_all{sub} = y;
    u_all{sub} = u;
    options_all{sub} = options;
end

%options for MFX
options_group.TolFun=1e-2;
options_group.MaxIter=50;
options_group.DisplayWin=0;
options_group.verbose=1;

priors_group.muPhi = 0; %temperature -- exp(phi(1))
priors_group.SigmaPhi = 10; %variance on temperature (before exponential transform)
priors_group.muTheta = zeros(dim.n_theta,1); %learning rate (alpha), selective maintenance (gamma) -- before logistic transform
priors_group.SigmaTheta =  1e1*eye(dim.n_theta); %variance of 10 on alpha and gamma
priors_group.muX0 = zeros(so.nbasis*so.hidden_states,1); %have PE and decay as tag-along states
priors_group.SigmaX0 = zeros(so.nbasis*so.hidden_states, so.nbasis*so.hidden_states); %have PE and decay as tag-along states
priors_group.a_vX0 = repmat(Inf, [1, so.nbasis*so.hidden_states]); %use infinite precision prior on gamma for X0 to treat as fixed (a = Inf; b = 0)
priors_group.b_vX0 = repmat(0, [1, so.nbasis*so.hidden_states]);

[p_sub, o_sub, p_group, o_group] = VBA_MFX_parallel(y_all, u_all, so.evo_fname, @g_sceptic, dim, options_all, priors_group, options_group);

%[p_sub, o_sub, p_group, o_group] = VBA_MFX(y_all, u_all, @h_sceptic_fixed_decay_fmri, @g_sceptic, dim, options_all, priors_group, options_group);

if is_alex
    cd ~/'Box Sync'/skinner/projects_analyses/SCEPTIC/mfx_analyses/
    if ~exist(so.model,'dir')
        mkdir(so.model)
        cd(so.model)
    end
end

%too huge to save into one .mat file
save('vba_mfx_results_psub.mat', 'p_sub', '-v7.3');
save('vba_mfx_results_pgroup.mat', 'p_group', '-v7.3');
save('vba_mfx_results_ogroup.mat', 'o_group', '-v7.3');
save('vba_mfx_results_osub.mat', 'o_sub', '-v7.3');
save('vba_mfx_results_settings.mat', 'priors_group', 'options_group', 'y_all', 'u_all', 'ids', '-v7.3');

delete(poolobj);
