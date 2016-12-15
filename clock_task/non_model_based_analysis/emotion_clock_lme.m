
%% Test the PIT idea: 
%  (1) in the instrumental context, appetitive outcome (reward) lengthens the subsequent RT 
%  because the relative value of staying in the trial increases relative to opportunity cost of doing something else.
%  (2) On the other hand, happy faces facilitate responding, shortening the RT, and vice versa for fearful faces.

%  get the data
cd('/Users/localadmin/Google Drive/skinner/projects_analyses/SCEPTIC/model_free_analyses');
load goldl; % only learnable contingencies, excluding missed and first trials of each run
load gold;

%% plot basic PE effect
%% looking just at learnable contingencies
figure(1); clf;
subplot(3,2,1)
scatter(gold.pemaxlag(gold.emotion_num == 1 & strcmpi(cellstr(gold.rewFunc),'DEV')),gold.rt(gold.emotion_num == 1 & strcmpi(cellstr(gold.rewFunc),'DEV'))); axis([-80 120 0 4000]); title('Happy, DEV'); ylabel('RT'); xlabel('Prediction error');
subplot(3,2,2)
scatter(gold.pemaxlag(gold.emotion_num == 1 & strcmpi(cellstr(gold.rewFunc),'IEV')),gold.rt(gold.emotion_num == 1 & strcmpi(cellstr(gold.rewFunc),'IEV'))); axis([-80 120 0 4000]); title('Happy, IEV'); ylabel('RT'); xlabel('Prediction error');
subplot(3,2,3)
scatter(gold.pemaxlag(gold.emotion_num == 0 & strcmpi(cellstr(gold.rewFunc),'DEV')),gold.rt(gold.emotion_num == 0 & strcmpi(cellstr(gold.rewFunc),'DEV'))); axis([-80 120 0 4000]); title('Scrambled, DEV'); ylabel('RT'); xlabel('Prediction error');
subplot(3,2,4)
scatter(gold.pemaxlag(gold.emotion_num == 0 & strcmpi(cellstr(gold.rewFunc),'IEV')),gold.rt(gold.emotion_num == 0 & strcmpi(cellstr(gold.rewFunc),'IEV'))); axis([-80 120 0 4000]); title('Scrambled, IEV'); ylabel('RT'); xlabel('Prediction error');
subplot(3,2,5)
scatter(gold.pemaxlag(gold.emotion_num == -1 & strcmpi(cellstr(gold.rewFunc),'DEV')),gold.rt(gold.emotion_num == -1 & strcmpi(cellstr(gold.rewFunc),'DEV'))); axis([-80 120 0 4000]); title('Fear, DEV'); ylabel('RT'); xlabel('Prediction error');
subplot(3,2,6)
scatter(gold.pemaxlag(gold.emotion_num == -1 & strcmpi(cellstr(gold.rewFunc),'IEV')),gold.rt(gold.emotion_num == -1 & strcmpi(cellstr(gold.rewFunc),'IEV'))); axis([-80 120 0 4000]); title('Fear, IEV'); ylabel('RT'); xlabel('Prediction error');


%%  test all effects with only intercept as random at each of the lower levels
lme_group_emotion_reward_abspelag = fitlme(goldl,'rt ~ 1 + condition_num + emotion*rewardlag + abspelag + rtvmaxlag + rtlag + rtlag2 +rtlag3 + emotion + (1|subject) + (1|run) + (1|run:subject)', 'DummyVarCoding', 'effects')
anova(lme_group_emotion_reward_abspe);
compare(lme_group_emotion_reward_abspe,lme_group_emotion_reward_abspelag)

%% plot interactions

% make a table for estimating predicted effects
t = table;
t.trial = [1:6]';
t.condition_num = goldl.condition_num(1:6);
t.abspelag = mean(goldl.abspelag)*ones(6,1);
t.emotion = {'fear'; 'fear'; 'scram'; 'scram'; 'happy'; 'happy'};
t.rewardlag = [0; 1; 0; 1; 0; 1];
t.rtvmaxlag = mean(goldl.rtvmax(~isnan(goldl.rtvmax)))*ones(6,1);
t.rtlag = mean(goldl.rt(~isnan(goldl.rt)))*ones(6,1);
t.rtlag2 = t.rtlag;
t.rtlag3 = t.rtlag;
t.subject = nominal(ones(6,1));
t.run = (ones(6,1));

% get predicted response

[yhat, yhatCI] = predict(lme_group_emotion_reward_abspelag,t);

% plot effects
fear = strcmpi(cellstr(t.emotion),'fear');
happy = strcmpi(cellstr(t.emotion),'happy');
scram = strcmpi(cellstr(t.emotion),'scram');
figure(99); clf;
errorbar(t.rewardlag(fear),yhat(fear),yhatCI(fear,1)-yhat(fear),yhatCI(fear,2)-yhat(fear)); hold on; 
errorbar(t.rewardlag(happy)+.05,yhat(happy),yhatCI(happy,1)-yhat(happy),yhatCI(happy,2)-yhat(happy)); 
errorbar(t.rewardlag(scram)-.05,yhat(scram),yhatCI(scram,1)-yhat(scram),yhatCI(scram,2)-yhat(scram)); hold off; 
legend('fear', 'happy', 'scrambled', 'Location', 'southeast');
ylabel('predicted response time')
ax = gca;
ax.XTick = [0 1]; ax.XTickLabel = {'Omission'; 'Reward'};
title('Reward by emotion interaction on clock task response times (error bars = 95% CI)'
