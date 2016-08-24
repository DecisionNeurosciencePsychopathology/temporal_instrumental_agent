%Vet TC generative model against Michael Frank's TC Feedback Aug2016

%From Michael Frank email 19Aug2016
K = 1500;
lambda = 0.2;
epsilon = 3000;
alphaG = 1;
alphaN = .5;
rho = 2500; % (or 0) 
nu = .3;

params = [lambda, epsilon, alphaG, alphaN, K, nu, rho]; %in order from original TC (used in forward)

nsubjs = 70;
rtbounds = [0 5000]; %define valid choice RTs (since agent is choosing freely)

%Per discussion with MJF, allow truly random draws from distribution over simulations.
%This eliminates perfect repeatability across simulations, but since we are averaging and smoothing
%across 'subjects,' the results should be convergent.
userngseed = 0;
ntrials = 50;

%complete one run of each contingency (order permuted within replications loop)
conds = {'DEV', 'IEV', 'CEV', 'CEVR'};

allRTpred = NaN(nsubjs, length(conds), ntrials);
allRTsmooth = NaN(nsubjs, length(conds), ntrials);
costs = NaN(nsubjs, length(conds));
smoother = 'loess';
span = 10; 

for i = 1:nsubjs
   condperm = conds(randperm(length(conds)));
   %setup priors for initial run
   priors.V = 0;
   priors.Go = 0;
   priors.NoGo = 0;
   priors.AvgRT = 2000; %average RT
   priors.FirstRT = 2000; %first choice

   %Gaussian variation in mean RT
   %priors.AvgRT = normrnd(2000, 400, 1);
   
   %Gaussian variation in first choice
   %priors.FirstRT = normrnd(2000, 400, 1);
   
   for j = 1:length(condperm)
       [cost, RTpred, ret] = TC_Alg_forward(params, priors, condperm{j}, -1, ntrials, rtbounds); %-1 for rngseed indicates not to use repeatable rng
       RTsmooth = smooth(RTpred, span, smoother);
       allRTpred(i, find(ismember(conds, condperm{j})), :) = RTpred; %find undoes the permutation so that indices follow cond order above
       allRTsmooth(i, find(ismember(conds, condperm{j})), :) = RTsmooth;
       priors.V = ret.V(length(ret.V)); %carry forweard value to next run
   end
end

%average over subjects
allRTavg = squeeze(mean(allRTpred, 1));
allRTsmoothavg = squeeze(mean(allRTsmooth, 1));

%plot results
figure(1); plot(allRTavg', 'LineWidth',2);
legend(conds); xlabel('Trial'); ylabel('Average RT across replications'); title(sprintf('TC results with rho = %.0f', rho));

figure(2); plot(allRTsmoothavg', 'LineWidth',2);
legend(conds); xlabel('Trial'); ylabel('Average RT across replications'); title(sprintf('Smoothed TC results with rho = %.0f', rho));
%print('TC_sims_withoutrho','-dpng','-r200')