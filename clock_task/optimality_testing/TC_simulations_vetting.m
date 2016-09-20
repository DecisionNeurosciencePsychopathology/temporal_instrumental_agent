%Vet TC generative model against Michael Frank's TC Feedback Aug2016
addpath('../')
%From Michael Frank email 19Aug2016
% simname='franksugg';
% K = 1500;
% lambda = 0.2;
% epsilon = 3000;
% alphaG = 1;
% alphaN = .5;
% rho = 2500; % (or 0) 
% nu = .3;

%best cherry picked parset I can find
% simname='good';
%params based on average in supplement of '09 (but tweaked to look reasonable)
% K = 1057.5;
% lambda = 0.14;
% epsilon = 3233.9;
% alphaG = 0.6;
% alphaN = 0.4;
% rho = 5360; % (or 0) 
% nu = 0.8;

simname='frank09_rho10k';
K = 1500;
lambda = 0.2;
epsilon = 3000;
alphaG = 0.3;
alphaN = 0.3;
nu = 0.2;
rho = 10000;

% simname='rho10k';
% %params based on average in supplement of '09 (but tweaked to look reasonable)
% K = 1100;
% lambda = 0.15;
% epsilon = 3200;
% alphaG = 0.6;
% alphaN = 0.4;
% rho = 10000; % (or 0) 
% nu = 0.7;

%from optimality
%           95148       2.9614       4.8925       60.664       4105.3       7719.1
% simname='optimized';
% K = 60.664;
% lambda = 0.00069646;
% epsilon = 95148;
% alphaG = 2.9614;
% alphaN = 4.8925;
% rho = 60.664;
% nu = 4105.3;   

params = [lambda, epsilon, alphaG, alphaN, K, nu, rho]; %in order from original TC (used in forward)

nsubjs = 500;
rtbounds = [0 5000]; %define valid choice RTs (since agent is choosing freely)

%Per discussion with MJF, allow truly random draws from distribution over simulations.
%This eliminates perfect repeatability across simulations, but since we are averaging and smoothing
%across 'subjects,' the results should be convergent.
userngseed = 0;
ntrials = 500;

%complete one run of each contingency (order permuted within replications loop)
conds = {'DEV', 'IEV', 'CEV', 'CEVR'};

allRTpred = NaN(nsubjs, length(conds), ntrials);
allRTsmooth = NaN(nsubjs, length(conds), ntrials);
costs = NaN(nsubjs, length(conds));
smoother = 'moving';
span = 15; 

for i = 1:nsubjs
   condperm = conds(randperm(length(conds)));
   %setup priors for initial run
   priors.V = 0;
   priors.Go = 0;
   priors.NoGo = 0;
   priors.AvgRT = 2500; %average RT
   priors.FirstRT = 1600; %first choice

   %Gaussian variation in mean RT
   %priors.AvgRT = normrnd(1700, 400, 1);
   
   %Gaussian variation in first choice
   priors.FirstRT = normrnd(1600, 100, 1);
   
   for j = 1:length(condperm)
       [cost, RTpred, ret] = TC_Alg_forward(params, priors, condperm{j}, -1, ntrials, rtbounds); %-1 for rngseed indicates not to use repeatable rng
       RTsmooth = smooth(RTpred, span, smoother);
       allRTpred(i, find(ismember(conds, condperm{j})), :) = RTpred; %find undoes the permutation so that indices follow cond order above
       allRTsmooth(i, find(ismember(conds, condperm{j})), :) = RTsmooth;
       priors.V = ret.V(length(ret.V)); %carry forward value to next run
   end
end

%average over subjects
allRTavg = squeeze(mean(allRTpred, 1));
allRTsmoothavg = squeeze(mean(allRTsmooth, 1));

%double smooth
allRTsmoothgroup = NaN(size(allRTsmoothavg));
for i = 1:size(allRTsmoothgroup, 1)
    allRTsmoothgroup(i,:) = smooth(allRTsmoothavg(i,:), 'moving', 10);
end

%plot results
figure(1); plot(allRTavg', 'LineWidth',2);
legend(conds); xlabel('Trial'); ylabel('Average RT across replications'); title(sprintf('TC results with rho = %.0f', rho));

figure(2); plot(allRTsmoothavg', 'LineWidth',2);
legend(conds); xlabel('Trial'); ylabel('Average RT across replications'); title(sprintf('Smoothed TC results with rho = %.0f', rho));

figure(3); plot(allRTsmoothgroup', 'LineWidth',2);
legend(conds); xlabel('Trial'); ylabel('Group smooth average RT across replications'); title(sprintf('Smoothed TC results with rho = %.0f', rho));
%print('TC_sims_withrho','-dpng','-r200')

save(sprintf('tcsims_%s', simname), 'allRTsmoothavg', 'allRTavg', 'allRTsmoothgroup', 'params', 'priors', 'ntrials', 'rtbounds');