%loads in subjects' data and fits SCEPTIC models using VBA;

% close all;
clear;
%curpath = fileparts(mfilename('fullpath'));

%behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
%%
%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
[~, me] = system('whoami');
me = strtrim(me);
if strcmp(os(1:end-2),'PCWIN')
    behavfiles = glob('C:\kod\temporal_instrumental_agent\clock_task\subjects\*.csv');
else
    if strcmp(me,'Alex')==1
        behavfiles = glob('/Users/localadmin/code/clock_smoothoperator/clock_task/subjects/*.csv');
        results_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
        group_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc';
    elseif strcmp(me(1:6),'dombax')==1
        behavfiles = glob('/Users/dombax/temporal_instrumental_agent/clock_task/subjects/*.csv');
        results_dir = '/Volumes/bek/vba_results/';
        group_dir = '/Users/dombax/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc';
        addpath(genpath('/Users/dombax/temporal_instrumental_agent/clock_task/'));
        addpath(genpath('/Users/dombax/code/'));
    elseif strcmp(me, 'wilsonj3')==1
        addpath(genpath('/Users/wilsonj3/matlab/temporal_instrumental_agent/clock_task/'));
        behavfiles = glob('/Users/wilsonj3/matlab/temporal_instrumental_agent/clock_task/subjects/*.csv');
        results_dir = '/Volumes/bek/vba_results/autocorrelation/';
        group_dir = '/Users/wilsonj3/matlab/temporal_instrumental_agent/clock_task/vba_results/autocorrelation';

    else
        behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
        addpath(genpath('/Users/dombax/temporal_instrumental_agent/clock_task/'));
        addpath(genpath('/Users/dombax/code/'));
    end
end

%% chose models to fit
modelnames = {'fixed' 'fixed_decay'};
%modelnames = {'fixed' 'fixed_uv' 'fixed_decay' 'kalman_softmax' 'kalman_processnoise' 'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'};
%modelnames = {'fixed' 'kalman_softmax' 'kalman_processnoise' 'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'};
%modelnames = {'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'}; %Rerun uncertainty models using corrected sigma update
% modelnames = {'kalman_logistic'};
%% set parameters
nbasis = 24;
multinomial = 1;
multisession = 0;
fixed_params_across_runs = 1;
fit_propspread = 0;
n_steps = 40;

u_aversion = 1; % allow for uncertainty aversion in UV_sum
saveresults = 1; %don't save to prevent script from freezing on Thorndike

% get ID list
id = NaN(length(behavfiles),1);

%% main loop
% L = NaN(length(modelnames),length(behavfiles));
% parpool
grp = struct([]);

fit_single_model = 0;
if fit_single_model
    model = 'fixed_decay'; % will run to get value and prediction errors.
    %     p = ProgressBar(length(behavefiles));
    parfor sub = 1:length(behavfiles)
        str = behavfiles{sub};
        id(sub) = str2double(str(isstrprop(str,'digit')));
        fprintf('Fitting subject %d id: %d \r',sub, id(sub))
        [posterior,out] = clock_sceptic_vba(id(sub),model,nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps,u_aversion);
        L(sub) = out.F;
        %gamma_decay(sub) = posterior.muTheta(2);
        tau(sub) = posterior.muTheta(1); %For fixed_uv
%         value(:,:,sub) = out.suffStat.muX;
        %         p.progress;
    end
    %     p.stop;
    cd(group_dir);
    filename = sprintf('nonMULTI_SHIFTED_U_check_grp_only_%s%d_nbasis%d_nsteps%d_uaversion%d',model,nbasis,n_steps, u_aversion);
    save(filename);
else
    %     %% continue where I left off earlie_r
    %     load L_n=64
    for m=1:length(modelnames)
        model = char(modelnames(m));
        %         p = ProgressBar(length(behavfiles));
        for sub=1:length(behavfiles)
            str = behavfiles{sub};
            
            %Really should just use regex here...
            if strcmp(me, 'wilsonj3')==1
                tmp_id = str(isstrprop(str,'digit'));
                id(sub) = str2double(tmp_id(2:end));
            else
                id(sub) = str2double(str(isstrprop(str,'digit')));
            end

            fprintf('Fitting %s subject %d \r',model,sub)
            %             p.progress;
            [posterior,out] = clock_sceptic_vba(id(sub),model,nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread,n_steps,u_aversion,saveresults);
%             cd(results_dir);
%             parsave(sprintf('output_%d',id(sub)),posterior,out);
            L(m,sub) = out.F;
        end
        %         p.stop;
        
    end
    cd(group_dir);
    filename = sprintf('SHIFTED_U_grp_L_%d_nbasis%d_nsteps%d_uaversion_not_allModels_fixed_prop_spread_entropy_elig_reward_variant',nbasis,n_steps, u_aversion);
    save(filename);
end
