%% Mixed-effects analyses of sceptic RT data

% data read 
cd('/Users/localadmin/Google Drive/skinner/projects_analyses/SCEPTIC/model_free_analyses');
load gold;

%% proprocessing
% define categorical variables
% gold = scepticfmribehaviorNov2016;
gold.exclude = gold.rt >3999 | gold.trial ==1;
gold.exclude_unlearnable = gold.rt >3999 | gold.trial ==1 | strcmpi(cellstr(gold.rewFunc),'CEVR') | strcmpi(cellstr(gold.rewFunc),'CEV');


gold.subject = categorical(gold.subject);
gold.LunaID = categorical(gold.LunaID);
gold.omission = categorical(gold.omission);
gold.omissionlag = categorical(gold.omissionlag);
gold.rewFunc = categorical(gold.rewFunc);
gold.emotion = categorical(gold.emotion);

gold.logrt = log(gold.rt);

%maybe 'nominal' rather than 'categorical'?
gold.omissionlag = nominal(gold.omissionlag);

% one example of how we make lagged variables:
gold.rewardlag2 = [NaN; gold.rewardlag(1:end-1)];
gold.rewardlag3 = [NaN; NaN; gold.rewardlag(1:end-2)];

gold.pemaxlag = [NaN; gold.pemax(1:end-1)];
gold.abspe = abs(gold.pemax);
gold.abspelag = abs(gold.pemaxlag);
gold.emotion_num = 1.*strcmpi(cellstr(gold.emotion),'happy')-1.*strcmpi(cellstr(gold.emotion),'fear');


%% derive within-run rtumax to avoid the problem of confounding with unique run intercept
% ugly, ugly code
first = 1;
gold.rtumax_run = NaN(size(gold.rtumax));
gold.rtumax_run_mean = NaN(size(gold.rtumax));

for unique_run = 1:length(gold.rtumax)/50;
    run_mean = mean(gold.rtumax(first:first+50));
    gold.rtumax_run(first:first+50) = gold.rtumax(first:first+50) - run_mean;
    gold.rtumax_run_mean(first:first+50) = run_mean;
    
    if unique_run<608
    first = first + 50;
    end
    unique_run
end
cols = fields(gold);
cols = cols(1:end-1);
goldl = gold(~gold.exclude_unlearnable,cols);
goldl.condition_num = 1.*strcmpi(cellstr(goldl.rewFunc),'IEV');
goldl.condition_num = nominal(goldl.condition_num);

gold.rtumax_runlag = [NaN; gold.rtumax_run(1:end-1)];
gold.rtumax_run_meanlag = [NaN; gold.rtumax_run_mean(1:end-1)];

gold.rtvmax_run = NaN(size(gold.rtumax));
for unique_run = 1:length(gold.rtumax)/50;
    run_mean = mean(gold.rtvmax(first:first+50));
    gold.rtumax_run = gold.rtvmax - run_mean;
    if unique_run<608
    first = first + 50;
    end
    unique_run
end


%% sanity check #1: inspect fields
gold(1:2,:)

%% sanity check #2: inspect rt distributions by condition
% scatterhist(gold.trial,gold.rt,'group',gold.rewFunc);figure(gcf)
figure(1); clf; boxplot(gold.rt,gold.rewFunc);
% 
% 
% % fit basic LM:
% lm_group = fitlm(gold,'rt ~ 1 + trial + subject + trial:subject + rewFunc + rewFunc:trial');
% 
% % fit basic LME:
% lme_group = fitlme(gold,'rt ~ 1 + rewFunc + trial + rewFunc:trial + (1 + rewFunc|subject) ');

%% explain RT in terms of reinforcement:
%% only intercept as random effect
lme_group_oneway = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtvmaxlag:trial + rtumaxlag + rtumaxlag:trial + (1|subject)')

%% all "meta-learning parameters" as random effects:
lme_group = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtvmaxlag:trial + rtumaxlag + rtumaxlag:trial + (1|subject) + (1+rtvmaxlag| subject) + (1+trial|subject) + (1+rtvmaxlag:trial|subject) + (1+rtumaxlag |subject) + (1+rtumaxlag:trial |subject)')

%% try to improve this model by removing nonsensical terms:
lme_group = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtvmaxlag:trial + rtumaxlag + rtumaxlag:trial + (1|subject) + (1| subject:run) + (rtvmaxlag:trial-1|subject) + (rtumaxlag -1|subject) + (rtumaxlag:trial-1|subject)')


% test naked effect of rtumaxlag
% lme_group_simplicissimo = fitlme(gold,'rt ~ 1 + rtvmaxlag + rtlag + rtlag2 + rtlag3 + rtumaxlag + (1|subject + 1| subject:run)')
lme_group_simplicissimo = fitlme(gold,'rt ~ 1 + rtvmaxlag + rtlag +  rtumaxlag + (rtumax|subject) + (rtumax|run) + (rtumax|run:subject)')


% add entropy
lme_group1 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtvmaxlag:trial + rtumaxlag + entropyH + rtumaxlag:trial + (1|subject) + (1| subject:run) + (rtvmaxlag:trial-1|subject) + (rtumaxlag -1|subject) + (rtumaxlag:trial-1|subject)')

% add PE
lme_group2 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtvmaxlag:trial + rtumaxlag + entropyH + pemax + rtumaxlag:trial + (1|subject) + (1| subject:run) + (rtvmaxlag:trial-1|subject) + (rtumaxlag -1|subject) + (rtumaxlag:trial-1|subject)')


% add rt_lag2
lme_group3 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtvmaxlag:trial + rtumaxlag + entropyH + pemax + rtumaxlag:trial + (1|subject) + (1| subject:run) + (rtvmaxlag:trial-1|subject) + (rtumaxlag -1|subject) + (rtumaxlag:trial-1|subject)')

% add rt_lag3
lme_group3 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtlag3 + rtvmaxlag:trial + rtumaxlag + entropyH + pemax + rtumaxlag:trial + (1|subject) + (1| subject:run) + (rtvmaxlag:trial-1|subject) + (rtumaxlag -1|subject) + (rtumaxlag:trial-1|subject)')

% did not work: add emotion, add emotion*PE, add entropyH*trial

% replace PE with scorelag
lme_group4 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtlag3 + rtvmaxlag:trial + rtumaxlag + entropyH + scorelag + rtumaxlag:trial + (1|subject) + (1| subject:run) + (rtvmaxlag:trial-1|subject) + (rtumaxlag -1|subject) + (rtumaxlag:trial-1|subject)')


% simplify random effects
lme_group5 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtlag3 + rtvmaxlag:trial + rtumaxlag + entropyH + rewardlag + rtumaxlag:trial + (1|subject) + (rtvmaxlag + rtumaxlag + rtvmaxlag:trial + rtumaxlag:trial | subject:run) ')

% estimate random effects at subject, not run level
lme_group6 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtlag3 + rtvmaxlag:trial + rtumaxlag + entropyH + rewardlag + rtumaxlag:trial + (1|subject:run) + (rtvmaxlag + rtumaxlag + rtvmaxlag:trial + rtumaxlag:trial | subject) ')

% same as before, but without random interactions
lme_group7 = fitlme(gold,'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtlag3 + rtvmaxlag:trial + rtumaxlag + entropyH + rewardlag + rtumaxlag:trial + (1|subject:run) + (rtvmaxlag + rtumaxlag + entropyH + rtlag + rewardlag| subject) ')

% add distance from edge
lme_group8 = fitlme(gold,  'rt ~ 1 + rtvmaxlag + trial + rtlag + rtlag2 + rtlag3 + rtvmaxlag:trial + rtumax_run_meanlag + rtumax_runlag + entropyH + rewardlag + distfromedgelag + rtumaxlag:trial + (1|subject:run) + (rtvmaxlag + rtumax_runlag + entropyH + rtlag + rewardlag|subject) ', 'CheckHessian', true)

% simple, non-RL model (now with lagged rewards [and respective random effects])
lme_group_nonRL = fitlme(gold,...
    'rt ~ 1 + trial + rtlag + rtlag2 + rtlag3 + rewardlag + rewardlag2 + rewardlag3 + distfromedgelag + (1|subject:run) + (rtlag + rewardlag| subject) ', 'CheckHessian', true)

% now add Value -- 
lme_group_nonRL_w_V_only = fitlme(gold,'rt ~ 1 + rtvmax + trial + rtlag + rtlag2 + rtlag3 + rewardlag + rewardlag2 + rewardlag3 + distfromedgelag + (1|subject:run) + (rtlag + rewardlag + rtvmax| subject) ', 'CheckHessian', true)


%% focus on emotion and prediction errors
% first rewards
lme_group_emotion_reward = fitlme(gold,'rt ~ 1 + emotion*rewardlag + emotion*rtvmaxlag + rtlag + rtlag2 +rtlag3 +  (emotion*rewardlag|subject) + (1|run) + (1|run:subject)', 'DummyVarCoding', 'effects') %, 'Exclude', find(gold.exclude))

%% switch go analyzing data without first trials, missed trials, and unlearnable contingencies

lme_group_emotion_reward_abspe = fitlme(goldl,'rt ~ 1 + condition_num + emotion*rewardlag + abspe*rewardlag + rtvmaxlag + rtlag + rtlag2 +rtlag3 + emotion*abspe + (1|subject) + (1|run) + (1|run:subject)', 'DummyVarCoding', 'effects')



% lme_group_emotion_pe = fitlme(gold,'rt ~ 1 + emotion*pemaxlag + emotion*rtvmaxlag + rtlag + rtlag2 +rtlag3 +  (emotion*pemaxlag|subject) + (1|run) + (1|run:subject)', 'DummyVarCoding', 'effects', 'Exclude', find(gold.exclude))



%% compare models
compare(lme_group_emotion_reward, lme_group_emotion_reward_abspe)
