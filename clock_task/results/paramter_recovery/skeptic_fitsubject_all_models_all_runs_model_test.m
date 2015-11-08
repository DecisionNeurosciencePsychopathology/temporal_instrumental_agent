% %JW: rmsearch was only minimizing costs at run level not subject level so
% %we hacked a solution to minimize cost locally in this function, instead of
% %minimizing in the actual optimizing function rmsearch_skeptic.m

function [cost] = skeptic_fitsubject_all_models_all_runs_model_test(params, sub, test_data, rngseeds,ntrials, nbasis, ntimesteps, modelname, sigma_noise,agent,optmat)
global count
count  = count + 1;
cost = 0; %initialize
for run=1;

    rts_obs = test_data(1,:);
    % apparently some RTs > 4000ms are recorded in the data: round them down to
    % 4000ms
    %rts_obs(rts_obs>4000) = 4000;
    
    %Need rew_rew as well!
    rts_rew = test_data(2,:);
    
    %cost = cost + skeptic_fitsubject_all_models(params, rts_obs', rew_obs', rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname);
    if strcmpi(agent,'frank')
        priors.V=0; %initialize expected value for first trial to prior (possibly from previous run)
        priors.Go=0; %initialize Go for first trial
        priors.NoGo=0; %initialize NoGo for first trial
        rtbounds = [1 5000]; %Don't let the agent choose 0
        cost = franktc_forRMsearch(params,priors, rts_obs, rts_rew,optmat,rngseeds, ntrials,rtbounds);
    else
        cost = clock_sceptic_agent_forRMsearch(params, rts_obs, rts_rew,'IEV', rngseeds,ntrials, nbasis, ntimesteps, 0,sigma_noise,agent);
    end                                  
    %sprintf('\n',params);
end

%fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,behav_sub.id, cost, modelname)
fprintf('Run #:%d \nSubject ID: %s\nSubject total cost: %.4f\nModel: %s\n\n',count,sub, cost, modelname)







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