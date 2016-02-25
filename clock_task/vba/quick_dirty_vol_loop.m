%loads in subjects' data and fits SCEPTIC models using VBA;

% close all;
clear variables;
%curpath = fileparts(mfilename('fullpath'));

global no_gamma

%behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
%%
%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    behavfiles = glob('C:\kod\temporal_instrumental_agent\clock_task\subjects\*.csv');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'Alex')==1
        behavfiles = glob('/Users/localadmin/code/clock_smoothoperator/clock_task/subjects/*.csv');
        results_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
        group_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc';
    elseif strcmp(me(1:6),'dombax')==1
        behavfiles = glob('/Users/dombax/temporal_instrumental_agent/clock_task/subjects/*.csv');
        results_dir = '/Volumes/bek/vba_results/sigma_volatility_variants';
        group_dir = '/Users/dombax/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc';
    else
        behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
        addpath(genpath('/Users/dombax/temporal_instrumental_agent/clock_task/'));
        addpath(genpath('/Users/dombax/code/'));
    end
end

%% chose models to fit
modelnames = {'kalman_sigmavolatility', 'fixed_uv'};

%% set parameters
nbasis = 4;
multinomial = 1;
multisession = 1;
fixed_params_across_runs = 1;
fit_propspread = 1;
n_steps = 10;

u_aversion = 1; % allow for uncertainty aversion in UV_sum
saveresults = 1; %don't save to prevent script from freezing on Thorndike

% get ID list
id = NaN(length(behavfiles),1);

%% main loop
% L = NaN(length(modelnames),length(behavfiles));
% parpool
grp = struct([]);

fit_single_model = 1;
for m=1:length(modelnames)
    for j = 1:2 %This loop should take care ot turning gamma on and off
%         if j==1
%             no_gamma = 0;
%         else
%             no_gamma = 1;
%         end
        model = char(modelnames(m));
        parfor sub = 1:length(behavfiles)
            str = behavfiles{sub};
            id(sub) = str2double(str(isstrprop(str,'digit')));
            fprintf('Fitting subject %d \r',sub)
            [posterior,out] = clock_sceptic_vba(id(sub),model,nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread, n_steps,u_aversion,saveresults);
            L(sub) = out.F;
            tau(sub) = posterior.muTheta(1);
            %         value(:,:,sub) = out.suffStat.muX;
            %         p.progress;
        end
        %     p.stop;
        cd(group_dir);
        filename = sprintf('grp_only_%s%d_nbasis%d_nsteps',model,nbasis,n_steps);
        %filename = sprintf('grp_only_%s%d_nbasis%d_nsteps%d_no_gamma',model,nbasis,n_steps, no_gamma);
        save(filename);
    end
end

