% %JW: rmsearch was only minimizing costs at run level not subject level so
% %we hacked a solution to minimize cost locally in this function, instead of
% %minimizing in the actual optimizing function rmsearch_skeptic.m

function [cost] = skeptic_fitsubject_all_models_all_runs_model_test(params, sub, test_data, rngseeds,ntrials, nbasis, ntimesteps,agent,sigma_noise,optmat,options)
global count
count  = count + 1;
cost = 0; %initialize


if isstruct(test_data)
    behav_sub = test_data;
    runs=unique(behav_sub.data.run);
else
    runs=1;
end


for run=1:length(runs);
    
    
    %Might need to check on the dimmensions of rts and rews
    if isstruct(test_data)
        cond = unique(behav_sub.data.rewFunc(behav_sub.data.run==runs(run)));
        %fprintf('\r%s subject#%d\r',char(cond),sub);
        %block = sprintf('%s_%d',char(cond),run);
        rts_obs = behav_sub.data.rt(behav_sub.data.run==runs(run));
        rts_obs = round(rts_obs/10)';
        % apparently some RTs > 4000ms are recorded in the data: round them down to
        % 4000ms
        rts_obs(rts_obs>4000) = 4000;
        rews = behav_sub.data.score(behav_sub.data.run==runs(run))';
        optmat=cond{:}; %just set this to condition for frank and qlearning
        %sigma_noise=0; %reset sigma noise
        sigma_noise = repmat(var(rews), 1, nbasis);
    else
        rts_obs = test_data(1,:);
        rews = test_data(2,:);
        if ismember(agent,{'franktc', 'franktc_fixed'}) || strcmpi(agent,'qlearning') || strcmpi(agent,'sarsa')
            optmat.sample = zeros(1,ntimesteps); %rest the optmat sampling back to 0 for each run
        end
        cond='IEV';
    end
    
    %cost = cost + skeptic_fitsubject_all_models(params, rts_obs', rew_obs', rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname);
    if ismember(agent,{'franktc', 'franktc_fixed'})
        priors.V=0; %initialize expected value for first trial to prior (possibly from previous run)
        priors.Go=0; %initialize Go for first trial
        priors.NoGo=0; %initialize NoGo for first trial
        rtbounds = [1 5000]; %Don't let the agent choose 0
        %%%%rts_obs = rts_obs.*10; %I commented this out when refitting the
        %%%%pseudo 'subjects' wanted to make sure rts were on the same
        %%%%scale.
        %optmat.sample = zeros(1,ntimesteps); %rest the optmat sampling back to 0 for each run
        cost = cost + franktc_forRMsearch(params,priors, rts_obs, rews,optmat,rngseeds, ntrials,rtbounds); %rews not really used
    elseif strcmpi(agent,'qlearning') || strcmpi(agent,'sarsa')
        options.episodeCount = ntrials;
        %could fix upstream but its late...
        if ~isstruct(test_data)
            rts_obs = test_data(1:ntrials,:)';
            rews = test_data(ntrials+1:end,:)';
        end
        %optmat.sample = zeros(1,ntimesteps); %rest the optmat sampling back to 0 for each run
        options.agent = agent;
        cost = cost + ClockWalking_3D_discountedEv_forRMsearch_subj_fitting(options,rts_obs,rews,optmat,rngseeds,params); %rews not really used
    else
        cost = cost + clock_sceptic_agent_forRMsearch(params, rts_obs, rews,cond, rngseeds,ntrials, nbasis, ntimesteps, 0,sigma_noise,agent);
    end                                  
    %sprintf('\n',params);
end

%fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,behav_sub.id, cost, modelname)
fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,sub, cost, agent)




% function [cost] = skeptic_fitsubject_all_models_all_runs(params, ~, behav_sub, rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname, ~)
% global count
% count  = count + 1;
% cost = 0; %initialize
% runs=unique(behav_sub.data.run);
% for run=1:length(runs);
%     %cond = unique(behav_sub.data.rewFunc(behav_sub.data.run==runs(run)));
%     %fprintf('\r%s subject#%d\r',char(cond),sub);
%     %block = sprintf('%s_%d',char(cond),run);
%     rts_obs = behav_sub.data.rt(behav_sub.data.run==runs(run));
%     % apparently some RTs > 4000ms are recorded in the data: round them down to
%     % 4000ms
%     rts_obs(rts_obs>4000) = 4000;
%     rew_obs = behav_sub.data.score(behav_sub.data.run==runs(run));
%     
%     cost = cost + skeptic_fitsubject_all_models(params, rts_obs', rew_obs', rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname);
%     %sprintf('\n',params);
% end
% 
% fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,behav_sub.id, cost, modelname)





%Old parameter recovery code
% function [cost] = skeptic_fitsubject_all_models_all_runs_model_test(params, sub, test_data, rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname, ~,sigma_noise)
% global count
% count  = count + 1;
% cost = 0; %initialize
% for run=1;
% 
%     rts_obs = test_data(1,:).*10;
%     % apparently some RTs > 4000ms are recorded in the data: round them down to
%     % 4000ms
%     %rts_obs(rts_obs>4000) = 4000;
%     
%     %Need rew_rew as well!
%     rts_rew = test_data(2,:);
%     
%     %cost = cost + skeptic_fitsubject_all_models(params, rts_obs', rew_obs', rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname);
%     cost = skeptic_fitsubject_all_models_new_peSelect(params, rts_obs, rts_rew, rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname,1,sigma_noise);
%     %sprintf('\n',params);
% end
% 
% %fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,behav_sub.id, cost, modelname)
% fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,sub, cost, modelname)