function [tGrid, fx] = distribPlots
% plots functions of the distributions used in qtask02.
% (for use in the cog sci talk and beyond)
% details:
%   1. functions are computed analytically from distrib params
%       -> not using empirically sampled delays from the experiment, as
%       in the cog sci paper
%       -> not using sampling methods
%   2. the two distributions may be plotted either on separate axes (better
%   for slides) or together on the same axes.
%   3. all functions need to account for truncation of the gp distribution

% basic parameters
tGrid = (0:0.01:20)'; % time grid (sec)
dist = {'unif', 'gp', 'beta'};
pars = {{0, 12}, {4, 5.75, 0}, {0.25, 0.25}};
scaling = [1, 1, 12];
trunc = [inf, 60, inf]; % upper bound at which the distribution is truncated
% for each distrib, G(t) = integral of t*f(t) 
% used for computing truncated mean (see comments below)
b = pars{1}{2};
Gt{1} = @(x) (x.^2)./(2*b);
k = pars{2}{1};
sigma = pars{2}{2};
Gt{2} = @(x) ((sigma+x)./(k-1)).*(1 + k.*x./sigma).^(-1/k);

% compute the functions of interest
fx = struct([]);
for d = 1:3 % for each distribution
    
    % find the cumulative probability corresponding to the point 
    % at which the distrib is truncated. the cdf at this point will
    % actually equal 1 due to truncation, and several functions will be
    % corrected accordingly.
    truncCDF = cdf(dist{d},trunc(d)./scaling(d),pars{d}{:});
    
    % probability density
    uncorrPDF = pdf(dist{d},tGrid./scaling(d),pars{d}{:})./scaling(d);
    fx(d).density = uncorrPDF./truncCDF; % correct for truncation
    fx(d).density(isinf(fx(d).density)) = 10; % deal with inf values
    
    % cumulative probability
    uncorrCDF = cdf(dist{d},tGrid./scaling(d),pars{d}{:});
    fx(d).cumul = uncorrCDF./truncCDF; % correct for truncation
    
    % hazard function
    fx(d).hazard = fx(d).density./(1-fx(d).cumul);

    % quartile boundaries
    quartBounds = zeros(1,4);
    for q = 1:4
        targProb = truncCDF*q/4;
        quartBounds(q) = scaling(d).*icdf(dist{d},targProb,pars{d}{:}); % desired percentile
    end
    fprintf('\t%s: quartile ceilings: %1.4f, %1.4f, %1.4f, and %1.4f sec\n',dist{d},quartBounds);
    
    % median-based prediction of remaining time
    %   (this will be shown in a dynamic plot)
    % first find the uncorrected CDF value corresponding to the
    % median-based prediction. this is the value halfway between the CDF at
    % the current time and the CDF at the truncation point
    predictionCDF = (uncorrCDF + truncCDF)./2;
    fx(d).medianAbs = scaling(d).*icdf(dist{d},predictionCDF,pars{d}{:}); % absolute prediction
    fx(d).median = fx(d).medianAbs - tGrid; % relative (i.e., time remaining)
    % the prediction is undefined for times above the distrib range
    rangeMax = trunc(d);
    if strcmp(dist{d},'unif'), rangeMax = min(rangeMax,pars{d}{2}); end
    if strcmp(dist{d},'beta'), rangeMax = scaling(d); end
    fx(d).median(tGrid>rangeMax) = nan;
    
    % expected return based on various quitting policies
    rwdLg = 0.05; % 20 points at 400 pts/$
    rwdSm = 0.0025; % 1 point at 400 pts/$
    iti = 0.8; % ITI duration
    totTime = 20*60; % 10 minutes to play
    % will vary with quitting policy:
    %   1. delay for the small reward (equals time waited before quitting)
    %   2. probability of large reward (equals the CDF)
    %   3. mean delay for large reward (more difficult to compute)
    delaySm = tGrid;
    
    %   pLg, probability of large reward, is the cdf of the reward
    %   distribution
    pLg = fx(d).cumul; 
    
    % SAMPLING approach to the bounded mean
    fprintf('sampling...');
    for t = 1:length(tGrid)
        tNow = tGrid(t);
        if mod(tNow,1)==0, fprintf('%d...',tNow); end
        if tNow==0
            delayLg(t,1) = 0;
        elseif tNow>rangeMax
            delayLg(t,1) = nan;
        else
            sampProb = uncorrCDF(t)*rand(1000,1); % uniform distrib of cumulative probabilities
            if strcmp(dist{d},'beta')
                sampProb(sampProb>0.999) = 0.999; % impose cutoff on beta dist for technical reasons
            end
            sampVal = scaling(d).*icdf(dist{d},sampProb,pars{d}{:});
            delayLg(t,1) = mean(sampVal);
        end
    end
    fprintf('\n');
    expRwd = pLg.*rwdLg + (1-pLg).*rwdSm;
    expDelay = pLg.*delayLg + (1-pLg).*delaySm + iti;
    sample(d).expReturn = totTime.*expRwd./expDelay;
    
    % use the sampling-based version for beta
    if strcmp(dist{d},'beta')
        fx(d).expReturn = sample(d).expReturn;
    else
        % analytical approach (described in comments below)
        delayLg = (Gt{d}(tGrid) - Gt{d}(0))./uncorrCDF;
        delayLg(tGrid>rangeMax) = nan;
        expRwd = pLg.*rwdLg + (1-pLg).*rwdSm;
        expDelay = pLg.*delayLg + (1-pLg).*delaySm + iti;
        fx(d).expReturn = totTime.*expRwd./expDelay;
    end
    
    % ANALYTICAL approach to the bounded mean
    % details - GP:
    %   the gp pdf f(x) = (1/sigma)*(1+k*(x/sigma))^(-1-1/k)
    %   (assuming location theta = 0)
    %   the integral of x*f(x) was computed using a calculator at
    %     http://integrals.wolfram.com
    %   it is G(x) = ((sigma+x)/(k-1))*(1 + k*x/sigma)^(-1/k)
    %   the mean truncated at T equals [G(T) - G(0)]/F(T)
    % details - UNIF:
    %   the mean truncated at T equals T/2
    %   this can be arrived at in the same way as above:
    %   G(x) = x^2/24 and F(x) = x/12
    %   [G(T) - G(0)]/F(T) = T/2
    % details - BETA:
    %   pdf f(x) = [x^(a-1) * (1-x)^(b-1)]/[normalization constant]
    %   we will ignore the normalization constant
    %   integral of f(x):
    %   F(x) = 
    
    % print info about the maximum earnings AND earnings at 12 sec
    [maxVal, maxIdx] = max(fx(d).expReturn);
    maxTime = tGrid(maxIdx);
    valAt12 = fx(d).expReturn(tGrid==12);
    fprintf('%s: max pay is $%1.2f at %1.3fsec; pay at 12sec is $%1.2f.\n',...
        dist{d},maxVal,maxTime,valAt12);
    
    
    
end

% create plots
fxNames = {'density', 'cumul', 'hazard', 'median', 'expReturn'};
fxYMax = [0.5, 1.1, 1, 22, 15];
fxYTick = {0:0.1:0.5, 0:0.2:1, 0:0.2:1, 0:5:20, 0:3:12};
fxYLabel = {'Probability density', 'Cumul. probability', 'Hazard rate', 'Time left (sec)', 'Return ($)'};
fxXLabel = [repmat({'Delay length (sec)'},1,4), {'Waiting policy (sec)'}];
col{1} = [0, 0, 1]; % unif condition: blue
col{2} = [1, 0, 0]; % gp condition: red 
col{3} = [0, 0.5, 0]; % beta condition: green 
for i = 1:length(fxNames)
    fname = fxNames{i};
    figure(i);
    fprintf('figure %d: %s\n',gcf,fname);
    for d = 1:3
        subplot(1,3,d);
        h = plot(tGrid,fx(d).(fname));
        % simulation overlay was included as a check (passed)
        % if strcmp(fname,'expReturn'), hold on; plot(tGrid,sample(d).expReturn,'r'); hold off; end
        set(gca,'XLim',[0, tGrid(end)],'YLim',[0, fxYMax(i)],...
            'XTick',0:5:tGrid(end),'YTick',fxYTick{i});
        set(gca,'Box','off','LineWidth',0.5,'FontSize',20);
        xlabel(fxXLabel{i});
        if d==1, ylabel(fxYLabel{i}); end
        set(h,'Color',col{d},'LineWidth',5);
        set(gcf,'Units','points','Position',[100, 230, 900, 275]); % 72 points per inch
        set(gca,'Units','points','Position',[25 + d*75 + (d-1)*200, 70, 200, 180]);
        
    end
end

% create plots for paper: multiple distributions on same axis. 
% goal is a figure 3 in wide, tiled 2x2 with plots: 1 in wide axes
fxNames = {'density', 'cumul', 'expReturn'};
fxYMax = [0.5, 1.1, 15];
fxYTick = {0:0.1:0.5, 0:0.2:1, 0:3:12};
fxYLabel = {'Probability density', 'Cumul. probability', 'Return ($)'};
fxXLabel = {'Delay length (sec)', 'Delay length (sec)', 'Waiting policy (sec)'};
col{1} = [0, 0, 0.5]; % unif condition: dark blue
col{2} = [1, 0.1, 0.1]; % gp condition: light red 
for i = 1:length(fxNames)
    fname = fxNames{i};
    figure(gcf+1);
    clf;
    hold on;
    for d = 1:3
        h = plot(tGrid,fx(d).(fname));
        set(gca,'XLim',[0, tGrid(end)],'YLim',[0, fxYMax(i)],...
            'XTick',0:5:tGrid(end),'YTick',fxYTick{i});
        set(gca,'Box','off','LineWidth',1,'FontSize',8);
        xlabel(fxXLabel{i});
        ylabel(fxYLabel{i});
        set(h,'Color',col{d},'LineWidth',2);
        set(gcf,'Units','points','Position',[100, 300, 144, 144]); % 72 points per inch
        axpos = get(gca,'Position');
        % set margins to 25%, axis dimensions to 50% of figure window
        set(gca,'Position',[0.25, 0.25, 0.5, 0.5]);
        
    end
end



