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

%modelnames = {'fixed' 'fixed_KL'};
%modelnames = {'v_process_noise_discounted'};
%modelnames = {'v_process_noise_discounted' 'fixed_rho'};
%modelnames = {'v_fixed_end_KL' 'v_fixed_rho_end_KL'};

%10/7/15 -- run this one first
%modelnames = {'fixed'};
%modelnames = {'fixed_rho'};
%modelnames = {'fixed_KL'};
 modelnames = {'v_processnoise'};
%modelnames = {'uv'};
%Then this one, make sure to change j!
%modelnames = {'v_fixed_rho_end_KL'};
%test_cases={'param_set1', 'param_set2', 'grw' 'param_set3'};
test_cases={'noise_params'};
ct=0;
%When fitting generated rts for specified paramters
%%%load('model_test_data2.mat')
%When computing std error of noise
load('temperature_test_data.mat')

%Set condition for seed
condition = 'IEV';

switch condition
    case 'IEV'
        cond_seed = 15;
        sigma_noise = temperature_test_data.(modelnames{:}).ret.set_1.sigma_nosie; %Should be the same for all data sets
end

%Set the beta
beta=0.1;
%beta=0.5;

for j = 7%:2:size(model_test_data,1)
    ct=ct+1;
% % %     test_data = model_test_data(j:j+1,:);
    test_case = test_cases{ct};
    %%%for sub = 1:2 %20 'subjects' all with same rts
    for sub = 1:length(fieldnames(temperature_test_data.(modelnames{:}).temperature_rts))
        
        %Grab the noisey generated rts and rews
        test_data(1,:)=temperature_test_data.(modelnames{:}).temperature_rts.(['set_' num2str(sub)]);
        test_data(2,:)=temperature_test_data.(modelnames{:}).temperature_rews.(['set_' num2str(sub)]);
        
        
        
        for i = 1:numel(modelnames)
            count=0;
            %% fit SKEPTIC RL model to each subject, using fixed 'optimal' parameters, with both
            % a full RT range and subject's actual (limited) RT range
            %         runs=unique(behav{sub}.data.run);
            %         for run=1:length(runs);
            
            if strcmpi(modelnames{i},'fixed')
                init_params_1 = 0.05;
                lower_bounds = 0;
                upper_bounds = 1;
            elseif strcmpi(modelnames{i},'fixed_KL')
                %Find alpha where 'fixed' did the best and set it to 'fixed discount'
                % initial param to hopefully cut down computation time
                %[~,idx] =min(fitted_vars.(modelnames{i-1}).(['subj_',sub]).cost_rmsearch);
                %init_params_1 = [fitted_vars.(modelnames{i-1}).(['subj_',sub]).fittedparameters_rmsearch(idx) 0 0]; %JW: 9/16/15 Changed kappa an lambda to 0, fittting problem...
                init_params_1 = [0.6 .007 .001];
                lower_bounds = [0 0 0];
                upper_bounds = [1 .3 .3]; %JW: 9/16/15 Changed from LB = -.1 UB = 1, bad fits...
            elseif strcmpi(modelnames{i},'v_process_noise_discounted') %JW: 9/17/15 added discounted process noise
                %Replace omega for alpha
                init_params_1 = [10 0 0];
                lower_bounds = [0 0 0];
                upper_bounds = [50 .125 .125];
                
            elseif strcmpi(modelnames{i},'v_processnoise') %JW: 9/17/15 added discounted process noise
                %Replace omega for alpha
                init_params_1 = 20; %Try making sigma noise a free parameter
                lower_bounds = 0;
                upper_bounds = 100;
                
            elseif strcmpi(modelnames{i},'fixed_rho') %JW: 9/17/15 added PE driven learning rates
                %Replace omega for alpha
                init_params_1 = [0.05 0.2];
                lower_bounds = [0 0];
                upper_bounds = [.5 .5];
                
            elseif strcmpi(modelnames{i},'v_fixed_end_KL') %JW: 9/29/15 added enduring models
                %Replace omega for alpha
                init_params_1 = [0.1 .1 .1];
                lower_bounds = [0 0 0];
                %upper_bounds = [.5 1 1];
                upper_bounds = [1 1 1];
                
            elseif strcmpi(modelnames{i},'v_fixed_rho_end_KL') %JW: 9/29/15 added enduring models
                %Replace omega for alpha
                init_params_1 = [0.1 0.3 .1 .1];
                lower_bounds = [0 0 0 0];
                %upper_bounds = [.5 .5 1 1];
                upper_bounds = [1 1 1 1];
            
            elseif strcmpi(modelnames{i},'uv') %JW: 10-14-15 uvSum test
                %Replace omega for alpha
                init_params_1 = 0.1;
                lower_bounds = 0;
                upper_bounds = 1;
            end
            
            %             cond = unique(behav{sub}.data.rewFunc(behav{sub}.data.run==runs(run)));
            %             fprintf('\r%s subject#%d\r',char(cond),sub);
            %             block = sprintf('%s_%d',char(cond),run);
            %             rts_obs = behav{sub}.data.rt(behav{sub}.data.run==runs(run));
            %             % apparently some RTs > 4000ms are recorded in the data: round them down to
            %             % 4000ms
            %             rts_obs(rts_obs>4000) = 4000;
            %             rew_obs = behav{sub}.data.score(behav{sub}.data.run==runs(run));
            % cond = behav{sub}.data.rewFunc(sub);
            %% fit the model with a full range of RTs
            % params = [.9165 .2261 .5]; %epsilon, prop_spread, spotlight
            %Need to add subject field
            fun = @(params) skeptic_fitsubject_all_models_all_runs_model_test([.2261 beta params],num2str(sub), test_data, [cond_seed 9 15 50], 24, 400, 0, 0, 500, modelnames{i},1,sigma_noise);
            %         optimoptions('fmincon')
            %         k = waitforbuttonpress;
            [fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).fittedparameters_rmsearch,fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).cost_rmsearch,...
                fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).exitflag_rmsearch,fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).xstart_rmsearch]=...
                rmsearch(fun, 'fmincon', init_params_1, lower_bounds, upper_bounds, 'initialsample', num_start_pts, 'options', opts,'plot','off');
            %Save the best cost and parameters associated with it.
            [best_cost,best_index] =min(fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).cost_rmsearch);
            fitted_vars.(test_case).(modelnames{i}).best_cost(sub,1) = best_cost;
            fitted_vars.(test_case).(modelnames{i}).best_parameters(sub,:) = fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).fittedparameters_rmsearch(best_index,:);
            
        end
        toc
    end
end
%Save file with timestamp
c = clock;
%JIC I forget to save the timestamp version
%save fitted_vars_newest_temp_model_test fitted_vars

%Fix timestamp save!
save(['fitted_vars_' c '.mat'],'fitted_vars', '-mat');

