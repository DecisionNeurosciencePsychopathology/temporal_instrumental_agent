
%% Test the PIT idea: (1) in the instrumental context, appetitive outcome (reward) triggers a longer subsequent RT 
lme_group_emotion_reward_abspe = fitlme(goldl,'rt ~ 1 + condition_num + emotion*rewardlag + abspe*rewardlag + rtvmaxlag + rtlag + rtlag2 +rtlag3 + emotion*abspe + (1|subject) + (1|run) + (1|run:subject)', 'DummyVarCoding', 'effects')
anova(lme_group_emotion_reward_abspe)
