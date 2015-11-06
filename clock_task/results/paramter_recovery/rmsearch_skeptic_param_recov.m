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
%modelnames = {'v_processnoise'};
%modelnames = {'uv'};
%Then this one, make sure to change j!
%modelnames = {'v_fixed_rho_end_KL'};
%modelnames = {'kalman_sigmavolatility'};
%modelnames = {'frank'};
% modelnames = {'kalman_uv_sum'};
% modelnames = {'kalman_uv_sum_kl'};
modelnames = {'fixedLR_softmax'};

%test_cases={'param_set1', 'param_set2', 'grw' 'param_set3'};
test_cases={'noise_params'};
ct=0;
%When fitting generated rts for specified paramters
%%%load('model_test_data2.mat')
%When computing std error of noise
load('param_recovery_test_data.mat')
load('param_recov.mat') %Load prev data
%Set condition for seed
condition = 'IEV';
rngseeds=[98 83 66 10];

%Set the beta
beta=0.1;
%beta=0.5;

%It might be more neat to use a struct like this, but it needs to be very
%adaptable, for instance if we want to use/not use propr spread and beta in
%the recovery, this will determine how many parameters we need and what the
%initial values of prop_spread and beta should be...
%a = initialize_agents_struct();

%lr_bounds=[.001 1];             %min/max values for all fixed learning rate parameters
prop_spread_bounds=[.001, .7];  %min/max values for prop_spread (gaussian temporal generalization)
beta_bounds=[.001, 2];
kl_bounds=[0 5];
lr_bounds=[.001 1]; 



%Start loop
for j = 1%:2:size(model_test_data,1)
    ct=ct+1;
% % %     test_data = model_test_data(j:j+1,:);
    test_case = test_cases{ct};
    fitted_vars.(test_case).(modelnames{1}).best_parameters = [];
    %%%for sub = 1:2 %20 'subjects' all with same rts
    for sub = 1:length(fieldnames(param_recovery_test_data.(modelnames{:}).ret))
        
        %Grab the noisey generated rts and rews
        %         test_data(1,:)=param_recovery_test_data.(modelnames{:}).ret.(['set_' num2str(sub)]).rts;
        %         test_data(2,:)=param_recovery_test_data.(modelnames{:}).ret.(['set_' num2str(sub)]).rew_i;
        test_data_1=param_recovery_test_data.(modelnames{:}).ret.(['set_' num2str(sub)]).rts;
        test_data_2=param_recovery_test_data.(modelnames{:}).ret.(['set_' num2str(sub)]).rew_i;
        test_data = [test_data_1; test_data_2];
        
        beta_init=.1;
        
%         init_params_1 = [0.2, 0.1]; %prop_spread, beta, tau (mix of U and V)
%         lower_bounds = [prop_spread_bounds(1), lr_bounds(1)];
%         upper_bounds = [prop_spread_bounds(2), lr_bounds(2)];
        
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
            elseif strcmpi(modelnames{i},'fixedLR_softmax')
                
                init_params_1 = [0.2, 0.1]; %prop_spread, beta, tau (mix of U and V)
                lower_bounds = [prop_spread_bounds(1), lr_bounds(1)];
                upper_bounds = [prop_spread_bounds(2), lr_bounds(2)];
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
                
            elseif strcmpi(modelnames{i},'kalman_sigmavolatility') 
                
                init_params_1 = [0.5, 0.8]; %phi, gamma prop spread beta are fixed!
                lower_bounds = [0, 0];
                upper_bounds = [10, 0.99];
                
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
                
            elseif strcmpi(modelnames{i},'frank')
                init_params_1 = [ 0.2, 3000, 0.3, 0.3, 1000, 0.1, 300 ]; %{'lambda', 'epsilon', 'alphaG', 'alphaN', 'K', 'nu', 'rho'};
                lower_bounds = [ 0, 0, 0.01, 0.01, 1, 0, 0 ];
                upper_bounds = [1, 100000, 5, 5, 5000, 5000, 10000 ];
            elseif strcmp(modelnames{i}, 'kalman_uv_sum');
                %Just fixing prop spread and beta
                prop_spread_init=0.62416337;           %tends to be pretty reasonable default for temporal generalization
                %beta_init=0.001524063;                    %default for temperature parameter.
                beta_init=0.1;                    %default for temperature parameter.
                %                 init_params_1 = 0.6; %prop_spread, beta, tau (mix of U and V)
                %                 lower_bounds = 0;
                %                 upper_bounds = 1;
                
                init_params_1 = [0.2, 0.6]; %prop_spread, beta, tau (mix of U and V)
                lower_bounds = [prop_spread_bounds(1), 0];
                upper_bounds = [prop_spread_bounds(2), 1];
            elseif strcmp(modelnames{i}, 'kalman_uv_sum_kl');
                %Just fixing prop spread and beta
                prop_spread_init=0.27298721;           %tends to be pretty reasonable default for temporal generalization
                %beta_init=1.2574299;                    %default for temperature parameter.
                beta_init=.1;
                %                 init_params_1 = 0.6; %prop_spread, beta, tau (mix of U and V)
                %                 lower_bounds = 0;
                %                 upper_bounds = 1;
                
                init_params_1 = [0.2, 0.6 0.1, 0.1]; %prop_spread, beta, tau (mix of U and V)
                lower_bounds = [prop_spread_bounds(1), 0, kl_bounds(1), kl_bounds(1)];
                upper_bounds = [prop_spread_bounds(2), 1, kl_bounds(2), kl_bounds(2)];
            end
            
            %% fit the model with a full range of RTs
            % params = [.9165 .2261 .5]; %epsilon, prop_spread, spotlight
            %Need to add subject field
            if strcmpi(modelnames{i},'frank')
                optmat=param_recovery_test_data.(modelnames{:}).ret.(['set_' num2str(sub)]).optmat;
                fun = @(params) skeptic_fitsubject_all_models_all_runs_model_test(params,num2str(sub), test_data, rngseeds,50,24,500,modelnames{i},0,modelnames{i},optmat);
            else
                sigma_noise_input = param_recovery_test_data.(modelnames{:}).ret.set_1.sigma_noise;
                %fun = @(params) skeptic_fitsubject_all_models_all_runs_model_test([prop_spread_init beta_init params],num2str(sub), test_data, rngseeds,50,24,500,modelnames{i},sigma_noise_input,modelnames{i});
                fun = @(params) skeptic_fitsubject_all_models_all_runs_model_test([params(1) beta_init params(2)],num2str(sub), test_data, rngseeds,50,24,500,modelnames{i},sigma_noise_input,modelnames{i});
            end

            %         optimoptions('fmincon')
            %         k = waitforbuttonpress;
            [fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).fittedparameters_rmsearch, fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).cost_rmsearch,...
                fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).exitflag_rmsearch, fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).xstart_rmsearch]=...
                rmsearch(fun, 'fmincon', init_params_1, lower_bounds, upper_bounds, 'initialsample', num_start_pts, 'options', opts,'plot','off');
            
            
            
            
            %Save the best cost and parameters associated with it.
            [best_cost,best_index] =min(fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).cost_rmsearch);
            
            
            fitted_vars.(test_case).(modelnames{i}).best_cost(sub,1) = best_cost;
            fitted_vars.(test_case).(modelnames{i}).best_parameters(sub,:) = fitted_vars.(test_case).(modelnames{i}).(['subj_' num2str(sub)]).fittedparameters_rmsearch(best_index,:);
            
            
            %Tried to use parfor loop on subject loop maybe come back to
            %it..
%             [fittedparameters_rmsearch, cost_rmsearch, exitflag_rmsearch, xstart_rmsearch]=...
%                 rmsearch(fun, 'fmincon', init_params_1, lower_bounds, upper_bounds, 'initialsample', num_start_pts, 'options', opts,'plot','off');
%             [best_cost,best_index] =min(cost_rmsearch);
%            best_cost(sub,1) = best_cost;
%            best_parameters(sub,:) = fittedparameters_rmsearch(best_index,:);
%            
%            temp_fitted_vars=save_to_struct(fitted_vars,modelnames{i},['subj_' num2str(sub)]);
           

        end
        toc
    end
end


%Save file with timestamp
c = clock;
%JIC I forget to save the timestamp version
%save fitted_vars_newest_temp_model_test fitted_vars


save param_recov fitted_vars 
save(['fitted_vars_' mat2str(c) '.mat'],'fitted_vars');

