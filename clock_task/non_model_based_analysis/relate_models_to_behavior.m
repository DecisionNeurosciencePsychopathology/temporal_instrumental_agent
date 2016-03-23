%% relate non-model-based results to relative model fits
% question 1: did the kalmans fit better on the 'optimal' subjects?
% question 2: did the kalmans fit worse on 'dislodged' subjects?

%% read in global max (6 consecutive trials) and dislodgement (6 consecutive off-max post global max) data
load('subject_rt_distribution_output');
subs = fields(subj_data);
for ct=1:length(subs)
    subj_data.global_max_reached(ct) = subj_data.(subs{ct}).global_max_reached;
    subj_data.dislodgement(ct) = subj_data.(subs{ct}).dislodgement_occurrence;
end

subj_data.dislodgement_binary = subj_data.dislodgement>0;

%% check distributions

close all;
figure(1); clf;
subplot(3,1,1);
hist(subj_data.global_max_reached,6)
subplot(3,1,2);
hist(subj_data.dislodgement,4)
subplot(3,1,3);
hist(subj_data.dislodgement_binary,2)

sum(subj_data.global_max_reached>3)

%% load fit data
load('models');
figure(2); clf; 
for ct=1:length(models.modelnames)
    subplot(length(models.modelnames),1,ct)
x = models.fixed_decay_adv_log_ratio(ct,:);
y = subj_data.global_max_reached;
    scatter(x,y)
    r = corr(x',y')
    

end

names = {'globs' 'disl' 'fix' 'fixUV' 'kSoft' 'kProc' 'kUV' 'kSigV' 'kLog' 'Q'};
% names(3:10) = models.modelnames([1:2 4:end]);

figure(3);clf;
corrplot([subj_data.global_max_reached' ...
    models.fixed_decay_adv_log_ratio([1:2 4:end],:)'],'VarNames',names([1,3:end]), 'type','Kendall','testR','on')
title('Global maxes reached vs. advantage of fixedDecay')

figure(4);clf;
corrplot([subj_data.dislodgement_binary(subj_data.global_max_reached>3)' ...
    models.fixed_decay_adv_log_ratio([1:2 4:end],subj_data.global_max_reached>3)'],'VarNames',names(2:end), 'type','Kendall','testR','on')
title('Dislodgement vs. advantage of fixedDecay')

