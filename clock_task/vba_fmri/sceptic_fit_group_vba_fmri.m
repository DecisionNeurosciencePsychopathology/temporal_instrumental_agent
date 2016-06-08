%loads in subjects' data and fits SCEPTIC models using VBA;

close all;
clear;
%curpath = fileparts(mfilename('fullpath'));

addpath(genpath('/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task'));
addpath('/storage/group/mnh5174_collab/temporal_instrumental_agent'); %has glob.m
addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox'));

behavfiles = glob('/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/subjects/*.csv');
group_dir = '/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/subjects';

%% set parameters
nbasis = 16; %16 basis functions total
multinomial = 1;
multisession = 0; %tend to get better fits for non-multisession
fixed_params_across_runs = 1;
fit_propspread = 0;
n_steps = 40;

saveresults = 1; %don't save to prevent script from freezing on Thorndike

% get ID list
id = NaN(length(behavfiles),1);

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

model = 'fixed_decay'; % will run to get value and prediction errors.
%     p = ProgressBar(length(behavefiles));
parfor sub = 1:length(behavfiles)
%for sub = 1:length(behavfiles)
    [~, str] = fileparts(behavfiles{sub});
    id(sub) = str2double(str(isstrprop(str,'digit')));
    fprintf('Fitting subject %d id: %d \r',sub, id(sub))
    [posterior,out] = clock_sceptic_vba_fmri(id(sub),model,nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps);
    L(sub) = out.F;
    %gamma_decay(sub) = posterior.muTheta(2);
    %tau(sub) = posterior.muTheta(1); %For fixed_uv
end
%     p.stop;
%cd(group_dir);
%filename = sprintf('nonMULTI_SHIFTED_U_check_grp_only_%s%d_nbasis%d_nsteps%d_uaversion%d',model,nbasis,n_steps,L);
%save(filename);

delete(poolobj);
