function [] = findhtdist()
% characterize probability distributions that will be used in experiment 2,
% and evaluate parameter values


globalScaling = 1; % may use this to slow or speed the entire design
tMax = 20*globalScaling; % in seconds
    % this is the length of the time interval within which conditions will
    % be compared. The uniform distribution extends to this length; the gp
    % distribution will contain delays that exceed this length. 
exptLength = 30*60;

whichDist = 'beta';
% whichDist = 'gp';
% whichDist = 'unif';

switch whichDist
    
    case 'beta'
        %%% beta distribution
        distpars = {0.3, 0.6};
        dist = 'beta';
        scaleFactor = tMax; % to scale the entire distribution
        cutoff = inf; % no cutoff
        
    case 'gp'
        %%% gp distribution
        distpars = {4, 5.75*3, 0};
        dist = 'gp';
        scaleFactor = 1*globalScaling;
        cutoff = 60*3;
        
    case 'unif'
        %%% uniform distribution
        distpars = {0, tMax};
        dist = 'unif';
        scaleFactor = 1;
        cutoff = inf; % no cutoff
        
end

% weibull distribution
% params = [4.8, .25]; % initial values for lambda, k
% dist = 'wbl';

% fixed parameters for timing and reward
consolation = 1/20; % proportion of reward obtained by quitting
iti = .8*globalScaling; % expressed in seconds
itiMin = tMax*consolation/(1 - 1.5*consolation);
    % with this value for ITI, a "maximum patience" policy is twice as
    % rewarding as a "maximum impatience" policy.
% if iti<itiMin, error('ITI is not long enough'); end

fprintf('distribution: %s( ',dist);
fprintf('%1.3f ',distpars{:});
fprintf(')\n');
fprintf('iti = %1.3f sec, consolation-prize ratio = %1.3f\n',iti,consolation);

% set the time grid
tvals = 0:0.1:tMax;
if tvals(end)~=tMax, error('rounding error creating time grid'); end

% proportion of the distribution that falls below the artificial cutoff
cdfCutoff = cdf(dist,cutoff./scaleFactor,distpars{:}); 

% compute expected time remaining as a function of time elapsed
% etr = calcETR(tvals,dist,distpars,cdfCutoff,scaleFactor);

% obtain distribution functions
cd = cdf(dist,tvals./scaleFactor,distpars{:});
pd = pdf(dist,tvals./scaleFactor,distpars{:});
cd_scaled = cd./cdfCutoff;
pd_scaled = pd./cdfCutoff;
haz = pd_scaled./(1 - cd_scaled);

% evaluate each policy i, viewed from the beginning of the interval
policyPayoff = nan(size(tvals));
policyCost = nan(size(tvals));
nSamp = 1000;
for i = 1:length(tvals)
    
    % print progress
    if mod(tvals(i),1)==0, fprintf('%ds...',tvals(i)); end
    
    % success rate under policy i
    p = cd_scaled(i); % probability a trial is rewarded
    policyPayoff(i) = p*1 + (1-p)*consolation; % average reward per trial
    
    % mean delay IF a trial is rewarded before the deadline
    range_p = cd(i)*rand(nSamp,1); % here we do NOT want the scaled version of the cdf
    s = warning('off','stats:betainv:NoConvergence');
            % icdf will sometimes produce a warning for high probabilities
            % with the beta distribution. This is ok, since it still
            % returns a 1 as it should.
    rangeSamps = icdf(dist,range_p,distpars{:});
    warning(s); % restore warning state
    meanIfRewarded = scaleFactor*mean(rangeSamps);
    % a trial's cost is measured in terms of its expected time in seconds
    % (a trial will either be rewarded or will last the full length)
    policyCost(i) = p*meanIfRewarded + (1-p)*tvals(i) + iti; % average per trial
end
fprintf('\n');

% summarize the benefit/cost ratio for all policies
payRatio = policyPayoff./policyCost;

% parameter-calibration results
defaultRate = 1/(0.5*tMax + iti);
fprintf('default pay ratio = %1.2f\n',defaultRate);
[peak peakInd] = max(payRatio);
totRewards = peak*exptLength;
fprintf('peak pay ratio = %1.2f (%1.2f times the default, %1.2f rewards total), reached at time %1.3f\n',...
    peak,peak/defaultRate,totRewards,tvals(peakInd));
comparWindow = tvals>0.1 & tvals<(tMax/2);
comparVals = payRatio - defaultRate;
comparVals(~comparWindow) = inf;
[matchVal matchInd] = min(abs(comparVals));
fprintf('default rate is best matched at time %1.3f (difference magnitude = %1.3f)\n',tvals(matchInd),matchVal);

% plot results
figure(1);
subplot(1,4,1);
% plot(tvals,etr,'k-','LineWidth',2);
title('expected time remaining');
subplot(1,4,2);
plot(tvals,pd_scaled,'k-','LineWidth',2);
title('pdf');
% set(gca,'YLim',[0, 1]);
subplot(1,4,3);
plot(tvals,cd_scaled,'k-','LineWidth',2);
title('cdf');
set(gca,'YLim',[0, 1]);
subplot(1,4,4);
plot(tvals(1:end-3),haz(1:end-3),'k-','LineWidth',2);
title('hazard');
% set(gca,'YLim',[0, 1]);

% plot rate-of-return results
figure(2); % costs and benefits
subplot(1,3,1);
plot(tvals,policyPayoff,'k-','LineWidth',2);
set(gca,'YLim',[0 1]);
title('policy payment per trial');
subplot(1,3,2);
plot(tvals,policyCost,'k-','LineWidth',2);
title('policy cost per trial (sec)');
subplot(1,3,3);
plot(tvals,payRatio,'k-','LineWidth',2);
title('benefit/cost ratio (rewards per sec)');
ylim = get(gca,'YLim');
set(gca,'YLim',[0, ylim(2)]);

end

% function to obtain the expected time remaining (for all values in tvals)
function [etr] = calcETR(tvals,dist,distpars,cdfCutoff,scaleFactor)
    
    % obtain cdf
    cd = cdf(dist,tvals./scaleFactor,distpars{:});
    
    % warn if the cdf does not begin at zero
    if cd(1)>0 
        disp('warning: distribution extends below zero');
    end

    % sample to estimate the mean time remaining at each point
    etr = zeros(size(tvals));
    nSamp = 10000; % number of samples for each mean
    for i = 1:length(tvals)
        cd_min = cd(i); % proportion of the component-1 distribution that has already gone by
        unifSamps = cd_min + (cdfCutoff-cd_min)*rand(nSamp,1); % uniform samples
        s = warning('off','stats:betainv:NoConvergence'); 
            % icdf will sometimes produce a warning for high probabilities
            % with the beta distribution. This is ok, since it still
            % returns a 1 as it should.
        distSamps = icdf(dist,unifSamps,distpars{:}); % convert to desired distribution
        warning(s); % restore warning state
        etr(i) = scaleFactor*mean(distSamps) - tvals(i); % remaining delay IF drawn from c1
    end

end













