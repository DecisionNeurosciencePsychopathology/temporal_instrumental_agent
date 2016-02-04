%loads in subjects' data and fits SCEPTIC models using VBA;

% close all;
clear variables;
%curpath = fileparts(mfilename('fullpath'));

%behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
%%
%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!

behavfiles = glob('../subjects/*.csv');
% results_dir = '/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results';
group_dir = 'E:\data\clock_task\vba\qlearning_vba_results';


%% chose models to fit
modelnames = {'qlearning'};

%% set parameters
nbasis = 4;
multinomial = 1;
multisession = 1;
fixed_params_across_runs = 1;
fit_propspread = 1;
n_steps = 40;
showfig = 0;
grp_flag = 1;


% get ID list
id = NaN(length(behavfiles),1);

%% main loop
% L = NaN(length(modelnames),length(behavfiles));
% parpool
grp = struct([]);

%Start parallel pool if not already running
try
    workers = 4;
    fprintf('\nStarting parallel pool with %d workers\n',workers);
    parpool('local', workers);
catch
    fprintf('\nParallel toolbox is most likely already runnning...\n');
end
    
fit_single_model = 1;
if fit_single_model
    model = 'qlearning'; %Define model
    parfor sub = 1:length(behavfiles)
        str = behavfiles{sub};
        id(sub) = str2double(str(isstrprop(str,'digit')));
        fprintf('Fitting subject %d \r',sub)
        [posterior,out] = clock_q_vba(id(sub),showfig, multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,grp_flag)
        L(sub) = out.F;
        tau(sub) = posterior.muTheta(1);
%         value(:,:,sub) = out.suffStat.muX;
        %         p.progress;
    end
    %     p.stop;
    cd(group_dir);
    filename = sprintf('grp_%s_nsteps_%d',model,n_steps);
    save(filename);
end
