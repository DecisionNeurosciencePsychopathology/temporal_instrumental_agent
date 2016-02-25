%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call Amoeba with our data
clear;
clc;
global count

options(1)=1;        % To display intermediate results use 1, otherwise use 0 (default)
options(2)=1e-3;     % Relative x-tolerance
options(3)=1e-3;     % Relative f-tolerance
options(14)=100;    % Max. number of f-evaluations per internal
fargs={};

tic
fitted_vars = struct;
% [fittedparameters_1,options]=simps('clock_logistic_operator',[.05 .95 1],[1 2 3],[options],[0.001 0.9 0.01],[0.2 .99 10],fargs{:});
% alpha_1=fittedparameters_1(1); lambda_1=fittedparameters_1(2); epsilon_1=fittedparameters_1(3);
% [cost_1, constr,value_all_1,value_hist_1]=clock_smooth_action_model(fittedparameters_1);

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    behavfiles = glob('C:\kod\temporal_instrumental_agent\clock_task\subjects\*.csv');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'Alex')==1
        behavfiles = glob('/Users/localadmin/clock_smoothoperator/clock_task/subjects/*.csv');
    else
        behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
    end
end



%epsilon is now scaled as the indifference point along mean(u_all) where agent switches from exploration to exploitation
%so should be a number between 0 and -1, roughly

%To look at fmincon options run this
%options = optimoptions('fmincon')

%fmincon_options = optimset(@fmincon);
%fmincon_options = optimoptions(@fmincon, 'UseParallel',true, 'Algorithm', 'active-set', 'DiffMinChange', 0.01);
%fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set', 'DiffMinChange', 0.005);
%JW 9/16/15 added maxFunEvals to quicken up the optimization
%options = optimoptions('fmincon');
options = optimoptions(@fmincon, 'UseParallel',true, 'Algorithm', 'active-set','Display','iter','DiffMinChange', 0.005);

opts = optimset('fmincon');
opts.LargeScale = 'off';
opts.Algorithm = 'active-set';
opts.Display = 'none';
%opts.MaxFunEvals = 1500;
%opts.TolFun = .001;
opts.DiffMinChange = .01;
%opts.MaxIter = 150;


num_start_pts=25;


%Choose which models to run
% x measns it was tested
% 1) 'fixedLR_softmax' x
% 2) 'fixedLR_egreedy' x
% 3) 'fixedLR_egreedy_grw'
% 4) 'asymfixedLR_softmax'
% 5) 'kalman_softmax' x
% 6) 'kalman_processnoise' x
% 7) 'kalman_sigmavolatility' x
% 8) 'kalman_uv_logistic' x
% 9) 'kalman_uv_sum' x
% 10) 'fixedLR_kl_softmax'
% 11) 'kalman_kl_softmax'
% 12) 'kalman_processnoise_kl'
% 13) 'kalman_uv_sum_kl' x
% 14) 'qlearning' x
% 15) 'sarsa'
% 16) 'franktc' x
% 17 'fixed_uv'
% 18 'frank_fixed'
a = initialize_stability_struct;

%Input model numbers from list
models = [17 18];


% a(models(i)).name = {};
% for i = 1:length(models)
%     a(models(i)).name = a(models(i)).name;
% end



%test_cases={'param_set1', 'param_set2', 'grw' 'param_set3'};
test_cases={'pseudo_subj_data_fitting' 'subj_fitting'};
ct=0;
%When fitting generated rts for specified paramters
%%%load('model_test_data2.mat')
%When computing std error of noise
load('param_recovery_test_data.mat')
load('param_recov.mat') %Load prev data
%Set condition for seed
condition = 'IEV';
rngseeds=[98 83 66 10];





%Start loop
ct=ct+1;

%Add another test case called subj_fitting
test_case = test_cases{1};


%Reset the best_parameters
for i = 1:length(models)
    fitted_vars.(test_case).(a(models(i)).name).best_parameters = [];
    fitted_vars.(test_case).(a(models(i)).name).best_cost = [];
end

if strcmpi(test_case,'subj_fitting')
    data = behavfiles;
elseif strcmpi(test_case,'pseudo_subj_data_fitting')
    data=fieldnames(param_recovery_test_data.(a(models(i)).name).ret);
end



%If anysubject breaks the fitting process
bad_apples = cell(length(data),1);

%CHECK ON SUB INITIAL NUMBER!!
for sub = 1:length(data) %Start with 8 since it died
    
    if strcmpi(test_case,'subj_fitting')
        %This guy is a bad apple skip him for now, he breaks process noise
        if sub==24
            continue
        end
        
        behav{sub}.data = readtable(behavfiles{sub},'Delimiter',',','ReadVariableNames',true);
        % write id
        fname = behavfiles{sub};
        idchars = regexp(fname,'\d');
        behav{sub}.id = fname(idchars);
        id = behav{sub}.id;
    else 
        id=sub; %Set id to sub by default
    end
    
    
    for i = 1:length(models)
        count=0;
        %% fit SKEPTIC RL model to each subject, using fixed 'optimal' parameters, with both
        % a full RT range and subject's actual (limited) RT range
        %         runs=unique(behav{sub}.data.run);
        %         for run=1:length(runs);
        
%         j = models(i);
        
    %If subject fitting set data as behav struct
    if strcmpi(test_case,'subj_fitting')
        test_data = behav{sub};
    else
        test_data_1=param_recovery_test_data.(a(models(i)).name).ret.(['set_' num2str(sub)]).rts;
        test_data_2=param_recovery_test_data.(a(models(i)).name).ret.(['set_' num2str(sub)]).rew_i;
        test_data = [test_data_1; test_data_2];
    end
        
        
        %% fit the model with a full range of RTs
        if ismember(a(models(i)).name,{'franktc', 'franktc_fixed'}) || strcmpi(a(models(i)).name,'qlearning') || strcmpi(a(models(i)).name,'sarsa')
            optmat=param_recovery_test_data.(a(models(i)).name).ret.(['set_' num2str(sub)]).optmat;
            a = initialize_stability_struct(id,test_data,rngseeds,0,optmat);
        else
            sigma_noise_input = param_recovery_test_data.(a(models(i)).name).ret.set_1.sigma_noise;
            a = initialize_stability_struct(id,test_data,rngseeds,sigma_noise_input);
        end
        
        
        try
        [fitted_vars.(test_case).(a(models(i)).name).(['subj_' num2str(id)]).fittedparameters_rmsearch, fitted_vars.(test_case).(a(models(i)).name).(['subj_' num2str(id)]).cost_rmsearch,...
            fitted_vars.(test_case).(a(models(i)).name).(['subj_' num2str(id)]).exitflag_rmsearch, fitted_vars.(test_case).(a(models(i)).name).(['subj_' num2str(id)]).xstart_rmsearch]=...
            rmsearch(a(models(i)).fun, 'fmincon', a(models(i)).init_params, a(models(i)).lower_bounds, a(models(i)).upper_bounds, 'initialsample', num_start_pts, 'options', opts,'plot','off');
        catch
            fprintf('Subject %d on agent %s broke rmsearch, logging...\n',sub,a(models(i)).name)
            bad_apples{sub,1} = {sub,a(models(i)).name};
            save bad_apples bad_apples
        end
        
        
        
        
        %Save the best cost and parameters associated with it.
        [best_cost,best_index] =min(fitted_vars.(test_case).(a(models(i)).name).(['subj_' num2str(id)]).cost_rmsearch);
        
        
        fitted_vars.(test_case).(a(models(i)).name).best_cost(sub,1) = best_cost;
        fitted_vars.(test_case).(a(models(i)).name).best_parameters(sub,:) = fitted_vars.(test_case).(a(models(i)).name).(['subj_' num2str(id)]).fittedparameters_rmsearch(best_index,:);
        
        
        %Tried to use parfor loop on subject loop maybe come back to
        %it..
        %             [fittedparameters_rmsearch, cost_rmsearch, exitflag_rmsearch, xstart_rmsearch]=...
        %                 rmsearch(fun, 'fmincon', init_params_1, lower_bounds, upper_bounds, 'initialsample', num_start_pts, 'options', opts,'plot','off');
        %             [best_cost,best_index] =min(cost_rmsearch);
        %            best_cost(sub,1) = best_cost;
        %            best_parameters(sub,:) = fittedparameters_rmsearch(best_index,:);
        %
        %            temp_fitted_vars=save_to_struct(fitted_vars,a(models(i)).name,['subj_' num2str(sub)]);
        
        
    end
    toc
    %In case of crashes?
    save temp_param_recov fitted_vars
end



%Save file with timestamp
c = clock;
%JIC I forget to save the timestamp version
%save fitted_vars_newest_temp_model_test fitted_vars


save param_recov fitted_vars
save(['fitted_vars_' mat2str(c) '.mat'],'fitted_vars');

