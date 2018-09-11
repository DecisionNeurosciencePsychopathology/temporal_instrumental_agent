%loads in subjects' data and fits SCEPTIC models using VBA;

close all;
clear;
%curpath = fileparts(mfilename('fullpath'));

os = computer;
[~, me] = system('whoami');
me = strtrim(me);
is_alex= strcmp(me,'Alex')==1 || strcmp(me,'dombax')==1;


%note that this function looks for 'sceptic_dataset' and 'sceptic_model'
%as environment variables so that this script can be scaled easily for batch processing
so = sceptic_validate_options(); %initialize and validate sceptic fitting settings


%% set environment and define file locations
if is_alex
  sceptic_repo=sprintf('/Users/%s/code/temporal_instrumental_agent/clock_task',me);
  addpath(genpath('~/code/VBA-toolbox')); %setup VBA
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
elseif strcmpi(so.dataset,'specc')
  behavfiles = glob([sceptic_repo, '/subjects/SPECC/*.csv']);
end

%extract IDs for record keeping
[~,fnames]=cellfun(@fileparts, behavfiles, 'UniformOutput', false);
ids=cellfun(@(x) char(regexp(x,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match')), fnames, 'UniformOutput', false);

%puts a nested cell in each element
%ids=regexp(fnames,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

so.output_dir = [sceptic_repo, '/vba_fmri/vba_out/', so.dataset, '/mfx/', so.model];
if ~exist(so.output_dir, 'dir'), mkdir(so.output_dir); end

%% setup parallel parameters
if is_alex
  if strcmp(me,'dombax')==1
      ncpus = 10;
  fprintf('defaulting to 10 cpus on Thorndike \n');
  else
  ncpus=4;
  fprintf('defaulting to 4 cpus on old iMac \n');
  end
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

ns = length(behavfiles);
y_all = cell(ns, 1);
u_all = cell(ns, 1);
options_all = cell(ns, 1);

n_t=NaN(1,ns);

for sub = 1:ns
  fprintf('Loading subject %d id: %s \n', sub, ids{sub});
  [data, y, u] = sceptic_get_data(behavfiles{sub}, so);
  [options, dim] = sceptic_get_vba_options(data, so);
  n_t(sub) = dim.n_t; % allow for variation in number of trials across subjects
  
  % populate data structures for VBA_MFX
  y_all{sub} = y;
  u_all{sub} = u;
  options_all{sub} = options;
end

% add the n_t vector to dim before passing to MFX
dim.n_t=n_t;

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
priors_group.a_vX0 = Inf([1, so.nbasis*so.hidden_states]); %use infinite precision prior on gamma for X0 to treat as fixed (a = Inf; b = 0)
priors_group.b_vX0 = zeros([1, so.nbasis*so.hidden_states]);

[p_sub, o_sub, p_group, o_group] = VBA_MFX_parallel(y_all, u_all, so.evo_fname, so.obs_fname, dim, options_all, priors_group, options_group);
%[p_sub, o_sub, p_group, o_group] = VBA_MFX(y_all, u_all, @h_sceptic_fixed_decay_fmri, @g_sceptic, dim, options_all, priors_group, options_group);

delete(poolobj);

%populate subject ids into output.options.inF structure since these are added to subject statistics below
for s=1:length(o_sub)
  o_sub{s}.options.inF.id = ids{s};
  o_sub{s}.options.inG.id = ids{s};

  %populate ffx parameters from o_group structure to p_sub structure for extraction
  p_sub{s}.ffx = o_group.initVBA.p_sub{s};
  p_sub{s}.ffx = rmfield(p_sub{s}.ffx, {'SigmaX', 'iQx', 'iQy'}); %large matrices not needed for anything (they use lots of disk space)
end

%create a structure with just the barebones useful parameters from subjects
s_all = cellfun(@(p, o) extract_subject_statistics(p, o), p_sub, o_sub, 'UniformOutput', false);

%output MFX results
if is_alex
  so.output_dir = '~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses';
end

[group_global, group_trial_level] = extract_group_statistics(s_all, ...
  sprintf('%s/%s_%s_mfx_sceptic_global_statistics.csv', so.output_dir, so.dataset, so.model), ...
  sprintf('%s/%s_%s_mfx_sceptic_trial_outputs_by_timestep.csv', so.output_dir, so.dataset, so.model));

%save the basis as a csv file
vnames = cellfun(@(x) strcat('Time_', num2str(x)), num2cell(1:so.ntimesteps), 'UniformOutput', false);
basis_mat = array2table(o_sub{1}.options.inF.gaussmat, 'VariableNames', vnames);

writetable(basis_mat, sprintf('%s/%s_%s_mfx_sceptic_basis.csv', so.output_dir, so.dataset, so.model));

%save group outputs
save(sprintf('%s/group_fits_%s_%s', so.output_dir, so.model, so.dataset), 'ids', 'so', 's_all', 'group_global', 'group_trial_level');

%too huge to save into one .mat file
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_psub.mat'], 'p_sub', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_pgroup.mat'], 'p_group', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_ogroup.mat'], 'o_group', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_osub.mat'], 'o_sub', '-v7.3');
save([so.output_dir, '/', so.dataset, '_', so.model, '_vba_mfx_results_settings.mat'], 'priors_group', 'options_group', 'y_all', 'u_all', 'ids', '-v7.3');