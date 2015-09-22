%This script will diagnose and compare the v_process noise and fixed
%learning rate models. Plots are cost comparision and overall value fitting
%comparisions

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

%% chose models to run

%modelnames = {'value_softmax' 'uv' 'v_discounted'};
% modelnames = {'v_discounted'};
%modelnames = {'uv_discounted'};
%modelnames = {'v_processnoise'};
%modelnames = {'v_processnoise' 'fixed'};
%modelnames = {'fixed_KL'};
%modelnames = {'fixed' 'fixed_KL' 'v_process_noise_discounted'};

%% Choose file
%filename = 'vars_noise_fixed_fixedKL';
%filename = 'fitted_vars_9-16-15';
% filename = 'fitted_vars_9-17-15';
% filename = 'fitted_vars_noise_run_1';
filename = 'fitted_vars_noiseDiscount_alphaOmega';
if exist([filename '.mat'], 'file')
    optimal_params = load([filename '.mat']);
end

modelnames = fieldnames(optimal_params.fitted_vars);

%% settings
trialplots = 1;
fit_params = 0;
compare_cost = 0;
%%

%wow, matlab has a useful function for mixed data types!!

%% write struct array of behavioral data with ids
for sub = 1:length(behavfiles)
    % write data
    behav{sub}.data = readtable(behavfiles{sub},'Delimiter',',','ReadVariableNames',true);
    % write id
    fname = behavfiles{sub};
    idchars = regexp(fname,'\d');
    behav{sub}.id = fname(idchars);
    
    %% fit SKEPTIC RL model to each subject, using fixed 'optimal' parameters, with both
    % a full RT range and subject's actual (limited) RT range
    runs=unique(behav{sub}.data.run);
    for run=1:length(runs);
        close all;
        cond = unique(behav{sub}.data.rewFunc(behav{sub}.data.run==runs(run)));
        fprintf('\r%s subject#%d\r',char(cond),sub);
        block = sprintf('%s_%d',char(cond),run);
        rts_obs = behav{sub}.data.rt(behav{sub}.data.run==runs(run));
        % apparently some RTs > 4000ms are recorded in the data: round them down to
        % 4000ms
        rts_obs(rts_obs>4000) = 4000;
        rew_obs = behav{sub}.data.score(behav{sub}.data.run==runs(run));
        % cond = behav{sub}.data.rewFunc(sub);
        %% fit the model with a full range of RTs
        % params = [.9165 .2261 .5]; %epsilon, prop_spread, spotlight
        
        
        % if fit_params
        
        beta = 100;
        for modelnum = 1:length(modelnames)
            modelname = char(modelnames(modelnum));
            if strcmpi(modelname,'value_softmax')
                params = [.2261 100]; %prop_spread, beta (temperature)
            elseif strcmpi(modelname,'uv') %% don't like this one anymore, don't use for parameter fitting
                params = [.2261 .9165 100]; %prop_spread, tau (relative weight of value vs. uncertainty), temperature
            elseif strcmpi(modelname,'v_processnoise') %scale process noise (learning rate) by abs(PE).
                params = [.2261 beta 10]; %prop_spread, temperature, scaling of PE into process noise
                opt_params = [.2261 beta optimal_params.subj_opt_vars.v_processnoise(sub,:)];
            elseif strcmpi(modelname,'v_process_noise_discounted') %scale process noise (learning rate) by abs(PE).
                params = [0.2261 beta 10 0 0]; %prop_spread, temperature, scaling of PE into process noise
                opt_params = [.2261 beta optimal_params.fitted_vars.v_process_noise_discounted.best_parameters(sub,:)];
            elseif strcmpi(modelname,'v_discounted')
                params = [.2261 .01 100]; %prop_spread, beta (temperature), kappa (discount factor)
            elseif strcmpi(modelname,'uv_discounted')
                %         params = [.2261 .01 .005 .9165]; %prop_spread, beta (temperature), kappa (discount factor), tau (V vs. U
                %         weight)
                %% NB -- higher than 'optimal' tau here
                params = [.2261 .01 .005 .95]; %prop_spread, beta (temperature), kappa (discount factor), tau
            elseif strcmpi(modelname, 'fixed')
                params = [0.2261 beta 0.2]; %prop_spread, beta (temperature), alpha (learning rate)
                %opt_params = [.2261 beta optimal_params.subj_opt_vars.fixed(sub,:)];
                
                %RMsearch params
                opt_params = [.2261 beta optimal_params.fitted_vars.fixed.best_parameters(sub,:)];
                
            elseif strcmpi(modelname, 'fixed_KL')
                params = [0.2261 beta 0.2 .0025 .0025]; %prop_spread, beta (temperature), alpha (learning rate)
                %opt_params = [.2261 beta optimal_params.subj_opt_vars.fixed_KL(sub,:)];
                
                %RMsearch Params
                opt_params = [.2261 beta optimal_params.fitted_vars.fixed_KL.best_parameters(sub,:)];
                
            elseif strcmpi(modelname, 'fixed_rho')
                params = [0.2261 beta 0.2 .02]; %prop_spread, beta (temperature), alpha (learning rate)
                %opt_params = [.2261 beta optimal_params.subj_opt_vars.fixed_KL(sub,:)];
                
                %RMsearch Params
                opt_params = [.2261 beta optimal_params.fitted_vars.fixed_rho.best_parameters(sub,:)];
            end
            
            %clock_logistic_fitsubject(params, rts_obs', rew_obs');
            range = 'full';
            %Quality check of models
            if modelnum == length(modelnames)
                next_model = 60;
            else
                next_model=0;
            end
            [~, ret, ~] = skeptic_fitsubject_all_models(params, rts_obs', rew_obs', [10 9 15 50], 24, 400, 0, 25, 400, modelname, run + modelnum);
            [~, ret, ~] = skeptic_fitsubject_all_models(opt_params, rts_obs', rew_obs', [10 9 15 50], 24, 400, trialplots, 25, 400, modelname, run+10 + modelnum);
            
            
            %Cost comparision
            if compare_cost
                [cost, ret, ~] = skeptic_fitsubject_all_models(opt_params, rts_obs', rew_obs', [10 9 15 50], 24, 400, 0, 25, 400, modelname, run);
                cost_h(sub,run,modelnum)=cost;
            end
            
            %% write p_choice fit statistics
            %behav{sub}.data.(sprintf('%s_p_chosen',modelname))(behav{sub}.data.run==runs(run),1) = ret.p_chosen';
        end
    end
end

%Plot cost comparision graph here
if compare_cost
    figure(42); clf;
    l_width = 2;
    for i = 1:run
        %Make a pretty graph
        subplot(4,2,i)
        plot(1:length(behavfiles),cost_h(:,i,1),'Linewidth', l_width) %v_noise
        hold on
        plot(1:length(behavfiles),cost_h(:,i,2), 'r', 'Linewidth', l_width) %fixed
        title(['Run: ' num2str(i)]);
        xlabel('Subjects')
        ylabel('Costs')
        
        %Which did better numerically
        higher_cost = cost_h(:,i,1)>cost_h(:,i,2);
        if sum(higher_cost)>round(sub/2)
            fprintf('v_processednoise had more subjects with a higher cost for run %i: %i\n',i,sum(higher_cost));
        else
            fprintf('fixed had more subjects with a higher cost for run %i: %i\n',i,sum(~higher_cost));
        end
    end
    legend({'v processnoise','fixed alpha'},'FontSize',8,'FontWeight','bold','Location','best')
    
end