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
% filename = 'fitted_vars_noiseDiscount_alphaOmega';
% if exist([filename '.mat'], 'file')
%     optimal_params = load([filename '.mat']);
% end

load('param_recov')
load('clock_options');
options = clock_options;
models = fieldnames(fitted_vars.subj_fitting);
pick = [2]; % skip Q and Frank for now;
models = models(pick);

%% settings
trialplots = 0;
fit_params = 0;
compare_cost = 1;
%%

nbasis=24;
ntimesteps=500;
reversal=0;
rngseeds=[98 83 66 10];

%wow, matlab has a useful function for mixed data types!!
start_with = 1;
%% write struct array of behavioral data with ids
for sub = start_with:length(fitted_vars.subj_fitting.fixedLR_softmax.best_cost);
    % write data
    behav{sub}.data = readtable(behavfiles{sub},'Delimiter',',','ReadVariableNames',true);
    % write id
    fname = behavfiles{sub};
    idchars = regexp(fname,'\d');
    behav{sub}.id = fname(idchars);
    
    %% fit SKEPTIC RL model to each subject, using fixed 'optimal' parameters, with both
    % a full RT range and subject's actual (limited) RT range
    runs=unique(behav{sub}.data.run);
    temp_cost=0;
    for run=1:length(runs);
        close all;
        cond = unique(behav{sub}.data.rewFunc(behav{sub}.data.run==runs(run)));
        fprintf('\r%s subject#%d\r',char(cond),sub);
        block = sprintf('%s_%d',char(cond),run);
        rts_obs = behav{sub}.data.rt(behav{sub}.data.run==runs(run));
        % apparently some RTs > 4000ms are recorded in the data: round them down to
        % 4000ms
        rts_obs = round(rts_obs/10);
        rts_obs(rts_obs>4000) = 4000;
        rew_obs = behav{sub}.data.score(behav{sub}.data.run==runs(run));
        % cond = behav{sub}.data.rewFunc(sub);
        %% fit the model with a full range of RTs
        % params = [.9165 .2261 .5]; %epsilon, prop_spread, spotlight
        
        
        % if fit_params
        
        beta = 100;
        beta=.1;
        for modelnum = 1:length(models)
            agent = char(models(modelnum));
            ret2=0; %set to 0 initially
            [opt_params, search_params, names] = getFittedParams(fitted_vars,sub,agent);
            
            %clock_logistic_fitsubject(params, rts_obs', rew_obs');
            range = 'full';
            %Quality check of models
            if modelnum == length(models)
                next_model = 60;
            else
                next_model=0;
            end
            ntrials=length(rew_obs);
            sigma_noise_input = repmat(var(rew_obs), 1, nbasis);
            
            %Need to add qlearning and frank
            if strcmp(agent,'qlearning')
                options.episodeCount = ntrials;
                options.agent = agent;
                [cost,~, ~,~,~,ret1]=ClockWalking_3D_discountedEv_forRMsearch_subj_fitting(options,rts_obs',rew_obs,cond,rngseeds,search_params);
                if compare_cost
                    [cost,~, ~,~,~,ret2]=ClockWalking_3D_discountedEv_forRMsearch_subj_fitting(options,rts_obs',rew_obs,cond,rngseeds,opt_params);
                    cost_h_qlearn(sub,run,1)=cost;
                end
            elseif  strcmp(agent,'franktc')
                priors.V=0; %initialize expected value for first trial to prior (possibly from previous run)
                priors.Go=0; %initialize Go for first trial
                priors.NoGo=0; %initialize NoGo for first trial
                rtbounds = [1 5000]; %Don't let the agent choose 0
                rts_obs_frank = rts_obs.*10;
                [cost, RTpred, ret1]=franktc_forRMsearch(search_params, priors,rts_obs_frank', rew_obs',cond, rngseeds, ntrials, rtbounds);
                
                if compare_cost
                    [cost,~,ret2]=franktc_forRMsearch(opt_params, priors,rts_obs_frank', rew_obs',cond, rngseeds, ntrials, rtbounds);
                    cost_h_frank(sub,run,1)=cost;
                end
            else
                
                [~,ret1] = clock_sceptic_agent_forRMsearch(search_params, rts_obs', rew_obs',cond, rngseeds, ntrials, nbasis, ntimesteps, reversal,sigma_noise_input,agent);
                %temp_cost = temp_cost+ clock_sceptic_agent_forRMsearch(search_params, rts_obs', rew_obs',cond, rngseeds, ntrials, nbasis, ntimesteps, reversal,sigma_noise_input,agent);
                % [~, ret, ~] = skeptic_fitsubject_all_models(opt_params, rts_obs', rew_obs', [10 9 15 50], 24, 400, 0, 25, 400, agent);
                
                %Cost comparision
                if compare_cost
                    [cost, ret2] = clock_sceptic_agent_forRMsearch(opt_params, rts_obs', rew_obs',cond, rngseeds, ntrials, nbasis, ntimesteps, reversal,sigma_noise_input,agent);
                    cost_h(sub,run,modelnum)=cost;
                end
            end
            
            
            %Plot subjects trial-wise
            fit_all_subs_all_models(agent,sub,ret1,ret2);
            
            
            %Investigating tau
            if strcmp(agent,'kalman_uv_sum')
                tau_range = [.1 .5 .95];
                for o = 1:length(tau_range)
                    search_params(3) = tau_range(o);
                    [~,ret_uv] = clock_sceptic_agent_forRMsearch(search_params, rts_obs', rew_obs',cond, rngseeds, ntrials, nbasis, ntimesteps, reversal,sigma_noise_input,agent);
                    fit_all_subs_all_models(agent,sub,ret_uv,ret2);
                end
            end
            
            
            %% write p_choice fit statistics
            %behav{sub}.data.(sprintf('%s_p_chosen',modelname))(behav{sub}.data.run==runs(run),1) = ret.p_chosen';
        end
    end
    stop=0;
end

%Plot cost comparision graph here
% if compare_cost
%     figure(42); clf;
%     l_width = 2;
%     for i = 1:run
%         %Make a pretty graph
%         subplot(4,2,i)
%         plot(1:length(behavfiles),cost_h(:,i,1),'Linewidth', l_width) %v_noise
%         hold on
%         plot(1:length(behavfiles),cost_h(:,i,2), 'r', 'Linewidth', l_width) %fixed
%         title(['Run: ' num2str(i)]);
%         xlabel('Subjects')
%         ylabel('Costs')
%
%         %Which did better numerically
%         higher_cost = cost_h(:,i,1)>cost_h(:,i,2);
%         if sum(higher_cost)>round(sub/2)
%             fprintf('v_processednoise had more subjects with a higher cost for run %i: %i\n',i,sum(higher_cost));
%         else
%             fprintf('fixed had more subjects with a higher cost for run %i: %i\n',i,sum(~higher_cost));
%         end
%     end
%     legend({'v processnoise','fixed alpha'},'FontSize',8,'FontWeight','bold','Location','best')
%
% end