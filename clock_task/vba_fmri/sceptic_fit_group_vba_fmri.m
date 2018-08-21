%loads in subjects' data and fits SCEPTIC models using VBA;

close all;
clear;
%curpath = fileparts(mfilename('fullpath'));

os = computer;
[~, me] = system('whoami');
me = strtrim(me);
is_alex=strcmp(me,'Alex')==1;

%note that this function looks for 'sceptic_dataset' and 'sceptic_model'
%as environment variables so that this script can be scaled easily for batch processing
so = sceptic_validate_options(); %initialize and validate sceptic fitting settings

%% set environment and define file locations
if is_alex
  sceptic_repo='~/code/temporal_instrumental_agent/clock_task';
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
elseif strcmpit(so.dataset,'specc')
  behavfiles = glob([sceptic_repo, '/clock_task/subjects/SPECC/*.csv']);
end

%extract IDs for record keeping
[~,fnames]=cellfun(@fileparts, behavfiles, 'UniformOutput', false);
ids=cellfun(@(x) char(regexp(x,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match')), fnames, 'UniformOutput', false);

%puts a nested cell in each element
%ids=regexp(fnames,'(?<=MEG_|fMRIEmoClock_)[\d_]+(?=_tc|_concat)','match'); %use lookahead and lookbehind to make id more flexible (e.g., 128_1)

so.output_dir = [sceptic_repo, '/vba_fmri/vba_out/', so.dataset];
if ~exist(so.output_dir, 'dir'), mkdir(so.output_dir); end

%% main loop
L = NaN(1,length(behavfiles));

ncpus=getenv('matlab_cpus');
if strcmpi(ncpus, '')
    ncpus=40;
    fprintf('defaulting to 40 cpus because matlab_cpus not set\n');
else
    ncpus=str2double(ncpus);
end

poolobj=parpool('local',ncpus); %just use shared pool for now since it seems not to matter (no collisions)

% p = ProgressBar(length(behavfiles));

parfor sub = 1:length(behavfiles)
  fprintf('Fitting subject %d id: %s \r', sub, ids{sub});
  
  [posterior, out] = clock_sceptic_vba_fmri(behavfiles{sub}, so);
  
  L(sub) = out.F;

  subj_id=ids{sub};

  %parsave doesn't work in recent MATLAB versions...
  m=matfile(sprintf('%s/sceptic_fit_%s_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d', ...
		  so.output_dir, ids{sub}, so.model, so.multinomial, so.multisession, ...
		  so.fixed_params_across_runs, so.u_aversion), 'writable',true);
  m.posterior=posterior; m.out=out; m.subj_id=ids{sub};
  
  %parsave(sprintf('%s/sceptic_fit_%s_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d', ...
  %		  so.output_dir, ids{sub}, so.model, so.multinomial, so.multisession, ...
  %    so.fixed_params_across_runs, so.u_aversion), posterior, out);%, subj_id);
  
  % gamma_decay(sub) = posterior.muTheta(2);
  % tau(sub) = posterior.muTheta(1); %For fixed_uv
end

%save group outputs for now
save(sprintf('%s/group_fits_%s_%s', so.output_dir, so.model, so.dataset), 'ids', 'L', 'so');

% p.stop;

delete(poolobj);
