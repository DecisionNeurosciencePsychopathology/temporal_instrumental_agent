function [so, poolobj, behavfiles] = sceptic_setup_environment(so)
%this function sets up the paths for sceptic, VBA, and expected behavior files based on the user's environment
%it also handles setup of the parpool depending on the user

%os = computer;
[~, me] = system('whoami');
me = strtrim(me);
is_alex= strcmp(me,'Alex')==1 || strcmp(me,'dombax')==1;
is_jiazhouchen=strcmp(me,'jiazhouchen')==1;

if is_alex
  if  strcmp(me,'Alex')
    sceptic_repo=sprintf('/Users/localadmin/code/temporal_instrumental_agent/clock_task');
      boxdir = sprintf('/Users/%s/Box',me);
  else
    sceptic_repo=sprintf('/Users/%s/code/temporal_instrumental_agent/clock_task', me);
  end
  
  addpath(genpath('~/code/VBA-toolbox')); %setup VBA
elseif is_jiazhouchen
  sceptic_repo=sprintf('/Users/jiazhouchen/Documents/UPMC/RStation/temporal_instrumental_agent/clock_task/');
  boxdir = sprintf('/Users/%s/Box',me);
else
  %ICS-ACI setup
  sceptic_repo='/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task';
  addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox')); %setup VBA
end

addpath(sceptic_repo); %has glob.m and setup_rbf.m

% identify subject behavior files
% group_dir=[sceptic_repo, '/subjects']; %not used for anything at the moment...

if strcmpi(so.dataset,'mmclock_meg')
  behavfiles = glob([sceptic_repo, '/subjects/mmclock_meg/*.csv']);
elseif strcmpi(so.dataset,'mmclock_fmri')
  behavfiles = glob([sceptic_repo, '/subjects/mmclock_fmri/*.csv']);
elseif strcmpi(so.dataset,'specc')
  behavfiles = glob([sceptic_repo, '/subjects/SPECC/*.csv']);
elseif strcmpi(so.dataset,'explore')
    rootdir = sprintf(fullfile(boxdir,'skinner','data','eprime','clock_reversal'));
    behavfiles = glob([rootdir, '/*/*.mat']);
end

if is_alex
  if strcmp(me,'dombax')
    so.output_dir = '/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses';
  else
    so.output_dir = '~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses';
  end
elseif is_jiazhouchen
  so.output_dir = '/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses';
else
  so.output_dir = [sceptic_repo, '/vba_fmri/vba_out/', so.dataset, '/mfx/', so.model];
  if ~exist(so.output_dir, 'dir'), mkdir(so.output_dir); end
end

% setup parallel parameters
if is_alex
  if strcmp(me,'dombax')==1
    ncpus = 12;
    fprintf('defaulting to 12 cpus on Thorndike \n');
  else
    ncpus=10;
    fprintf('defaulting to 10 cpus on iMac Pro \n');
  end
  poolobj = gcp('nocreate');
  if isempty(poolobj)
    poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)
  end
elseif is_jiazhouchen
  ncpus=4;
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

end
