%loads in subjects' data and fits SCEPTIC models using VBA;

% close all;
clear;
%curpath = fileparts(mfilename('fullpath'));

behavfiles = glob('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/*.csv');
addpath(genpath('/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task'));

%% chose models to fit
%modelnames = {'fixed' 'fixed_uv' 'fixed_decay' 'kalman_softmax' 'kalman_processnoise' 'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'};
%modelnames = {'fixed' 'kalman_softmax' 'kalman_processnoise' 'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'};
modelnames = {'fixed_uv' 'kalman_uv_sum'}; %Rerun uncertainty models using corrected sigma update

%% set parameters
nbasis = 16; %16 basis functions total
multinomial = 1;
multisession = 1;
fixed_params_across_runs = 1;
fit_propspread = 1;
n_steps = 40;

saveresults = 1; %don't save to prevent script from freezing on Thorndike

% get ID list
id = NaN(length(behavfiles),1);

%% main loop
% L = NaN(length(modelnames),length(behavfiles));
% parpool
grp = struct([]);


model = 'fixed_decay'; % will run to get value and prediction errors.
%     p = ProgressBar(length(behavefiles));
parfor sub = 1:length(behavfiles)
%for sub = 1:length(behavfiles)
    str = behavfiles{sub};
    id(sub) = str2double(str(isstrprop(str,'digit')));
    fprintf('Fitting subject %d id: %d \r',sub, id(sub))
    [posterior,out] = clock_sceptic_vba_fmri(id(sub),model,nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps);
    L(sub) = out.F;
    %gamma_decay(sub) = posterior.muTheta(2);
    tau(sub) = posterior.muTheta(1); %For fixed_uv
end
%     p.stop;
cd(group_dir);
filename = sprintf('nonMULTI_SHIFTED_U_check_grp_only_%s%d_nbasis%d_nsteps%d_uaversion%d',model,nbasis,n_steps);
save(filename);
