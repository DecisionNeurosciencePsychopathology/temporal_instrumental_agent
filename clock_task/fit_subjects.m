%script to load in subjects' data and fit using the logistic operator
%behavfiles = dir( fullfile('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files', '*.csv'));
%curpath = fileparts(mfilename('fullpath'));

%behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
%%
[~, me] = system('whoami');
me = strtrim(me);
if strcmp(me,'Alex')==1
  behavfiles = glob('/Users/localadmin/clock_smoothoperator/clock_task/subjects/*.csv');
else    
  behavfiles = glob('/Users/michael/Data_Analysis/clock_analysis/fmri/behavior_files/*.csv');
end


%%
%read header
% fid = fopen(behavfiles(1).name, 'r');
% head = textscan(fid, '%s', 1);
% fclose(fid);
% m = csvread(behavfiles(1).name, 1, 0);

%% settings
trialplots = 1;
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


modelname = 'value_softmax';
if strcmpi(modelname,'value_softmax')
params = [.2261 .01]; %prop_spread, spotlight
else
end

%clock_logistic_fitsubject(params, rts_obs', rew_obs');
range = 'full';
[~, ret, ~] = skeptic_fitsubject_all_models(params, rts_obs', rew_obs', [10 9 15 50], 24, 400, trialplots, 25, 400, modelname);



%version with spotlight
%[~, ret, ~] = skeptic_fitsubject(params(1:3), rts_obs', rew_obs', [10 9 15 50], 24, 400, trialplots, cond);

% %% write predicted RTs -explore and -exploit
 behav{sub}.data.full_rtpred_explore(behav{sub}.data.run==runs(run),1) = ret.rts_pred_explore(1:50)';
 behav{sub}.data.full_rtpred_exploit(behav{sub}.data.run==runs(run),1) = ret.rts_pred_exploit(1:50)';


%% test fits with regression
fprintf('\r**************\r%s range\r',range);
behav{sub}.stats_full.(block)=test_reg(rts_obs,ret, modelname);
fprintf('**************\r');
% 
% %% Let's also try with a representational basis restricted to the subject's RT range
% range = 'limited';
% minrt = round(min(rts_obs(rts_obs>0))./10);
% maxrt = round(max(rts_obs)./10);
% [~, ret, ~] = skeptic_fitsubject_all_models(params, rts_obs', rew_obs', [10 9 15 50], 24, 400, trialplots, minrt, maxrt, modelname);
% % 
% %  behav{sub}.data.limited_rtpred_explore(behav{sub}.data.run==runs(run),1) = ret.rts_pred_explore(1:50)';
% %  behav{sub}.data.limited_rtpred_exploit(behav{sub}.data.run==runs(run),1) = ret.rts_pred_exploit(1:50)';
% 
% % %% write predicted RTs for representational range limited to subject's RT range
% % behav{sub}.data.lim_rtpred_explore(behav{sub}.data.run==runs(run)) = ret.rts_pred_explore(1:50)';
% % behav{sub}.data.lim_rtpred_exploit(behav{sub}.data.run==runs(run)) = ret.rts_pred_exploit(1:50)';
% 
% %% test with regression
% fprintf('\r**************\r%s range\r',range);
% behav{sub}.stats_lim.(block)=test_reg(rts_obs,ret, modelname);
% fprintf('**************\r');
end


end
save behav
%% let's see if R2 improves for limited vs. full range
%% more importantly, are betas for RT_pred(explore, exploit) different from 0?
for sub=1:length(behavfiles)
    conds = fields(behav{sub}.stats_full);
for run = 1:length(conds)
r2full(sub,run) = behav{sub}.stats_full.(conds{run}).stats(1);
r2lim_range(sub,run) = behav{sub}.stats_lim.(conds{run}).stats(1);

b_exploit_full(sub,run) = behav{sub}.stats_full.(conds{run}).b(1);

b_L_CI_exploit_full(sub,run) = behav{sub}.stats_full.(conds{run}).bint(1,1);
b_U_CI_exploit_full(sub,run) = behav{sub}.stats_full.(conds{run}).bint(1,2);

b_explore_full(sub,run) = behav{sub}.stats_full.(conds{run}).b(2);

b_L_CI_explore_full(sub,run) = behav{sub}.stats_full.(conds{run}).bint(2,1);
b_U_CI_explore_full(sub,run) = behav{sub}.stats_full.(conds{run}).bint(2,2);


b_exploit_lim(sub,run) = behav{sub}.stats_lim.(conds{run}).b(1);
b_explore_lim(sub,run) = behav{sub}.stats_lim.(conds{run}).b(2);
end
end

%% sanity check plots
% clf
% figure(99);
% errorbar(1:length(behavfiles),b_exploit_full(:,2),b_exploit_full(:,2)-b_L_CI_exploit_full(:,2),b_U_CI_exploit_full(:,2)-b_exploit_full(:,2));hold on;
% errorbar(1:length(behavfiles),b_explore_full(:,2),b_explore_full(:,2)-b_L_CI_explore_full(:,2),b_U_CI_explore_full(:,2)-b_explore_full(:,2),'r');hold off;
% 
% figure(98)
% scatter(b_exploit_full(:,2),b_explore_full(:,2));
%% testing range limitation

[h,p,CI,stats] = ttest(r2full,r2lim_range);
%no, R2 does not seem to improve

%% are betas different from 0?
% this should really be done using CIs stored above!!!
[h,p,CI,stats] = ttest(b_exploit_lim,0);
[h,p,CI,stats] = ttest(b_explore_lim,0);
% strangely, 
figure(1);clf;
subplot(4,1,1);
hist(b_explore_full,40);
xlabel('RTexploRE regression coefficients; color: run')
ylabel('Subject count')
title('Full RT range results');
subplot(4,1,2)
hist(b_exploit_full,200)
axis([-5 5 0 200])
xlabel('RTexploIT regression coefficients; color: run')
ylabel('Subject count')
title('Full RT range results WITHOUT OUTLIERS');

subplot(4,1,3);
hist(b_explore_full,40);
xlabel('RTexploRE regression coefficients; color: run')
ylabel('Subject count')
title('Limited RT range results');
subplot(4,1,4)
hist(b_exploit_full,200)
 axis([-5 5 0 200])
xlabel('RTexploIT regression coefficients; color: run')
ylabel('Subject count')
title('Limited RT range results WITHOUT OUTLIERS');


%% export data to R

for sub=1:length(behavfiles)
    writetable(behav{sub}.data,sprintf('%s.csv',behav{sub}.id));
end
%just example first run, single subject
%sub = 53;
%%
%function [cost,mov,ret] = skeptic_fitsubject(params, rt_obs, rew_obs, rngseeds, nbasis, ntimesteps, trial_plots, cond, minrt)

fmincon_options = optimoptions(@fmincon, 'UseParallel',false, 'Algorithm', 'active-set');%, 'DiffMinChange', 0.001);

init_params_1 = [.02 .9877 -.06];
lower_bounds = [0.001 0.9 -1];
upper_bounds = [0.2 .999 0];

[fittedparameters_fmincon, cost_fmincon, exitflag_fmincon] = fmincon(@(params) clock_logistic_fitsubject(params, rts_obs', rew_obs'), init_params_1, [], [], [], [], lower_bounds, upper_bounds, [], fmincon_options);

[~, fitted_object] = clock_logistic_fitsubject(-.09, rts_obs', rew_obs');
fitted_object.cost_total
[~, fitted_object] = clock_logistic_fitsubject(-.06, rts_obs', rew_obs');
fitted_object.cost_total

plot(1:50, fitted_object.rts_obs)
hold on;
plot(1:50, fitted_object.rts_pred, 'b')
