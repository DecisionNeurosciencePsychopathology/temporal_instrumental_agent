function [cost, RTpred, ret]=TC_Alg_forward(params, priors, cond, rngseeds, ntrials, rtbounds)
% Time clock R-L algorithm developed by Michael Frank
% adapted to be generative/forward agent only for testing
% only the beta distribution version here (not Kalman)
%
%
% inputs:
%    params:   vector of model parameters used to fit data
%    cond:  reward contingency, 1=CEV; 2=CEVR; 3=DEV; 4=IEV
%    ntrials:  number of trials to execute
%
%
% params is an 7 x 1 column vector of model parameters to be used to fit behavior.

if nargin < 4, rngseeds=31; end
if nargin < 5, ntrials=50; end
if nargin < 6, rtbounds=[0 5000]; end

%fprintf('Fitting cond %s\n', cond);

global rew_rng_state;

if rngseeds == -1
    userngseed = 0; %truly random draws
else
    %initialize state for a repeatable random number generators using different seeds (used by RewFunction)
    rew_rng_seed=rngseeds(1);
    
    rng(rew_rng_seed);
    rew_rng_state=rng;
    userngseed = 1;
end

%defaults for FirstRT and AvgRT if not set in priors
if ~isfield(priors, 'FirstRT'), priors.FirstRT = sum(rtbounds)/2; end %midpoint
if ~isfield(priors, 'AvgRT'), priors.AvgRT = sum(rtbounds)/2; end %midpoint

RTpred = NaN(ntrials, 1); %vector of predicted RTs
Reward = NaN(ntrials, 1); %vector of obtained rewards
V = NaN(ntrials, 1); %state-value function (expected value)
ev = NaN(ntrials, 1); %expected value of choices (for cost function)
Go = NaN(ntrials, 1);
NoGo = NaN(ntrials, 1);

usestruct=0;
if isstruct(cond)
    cstruct = cond;
    cond = cstruct.name;
    usestruct=1;
end

RTpred(1) = priors.FirstRT; %typically start in middle

%obtain first outcome
if usestruct
    [Reward(1), ev(1), cstruct] = getNextRew(round(RTpred(1)/10)+1, cstruct);
else
    [Reward(1), ev(1)] = RewFunction(RTpred(1), cond, userngseed);
end

V(1) = priors.V; %initialize expected value for first trial to prior (possibly from previous run)
Go(1) = priors.Go; %initialize Go for first trial
NoGo(1) = priors.NoGo; %initialize NoGo for first trial

%V_fast = V(1); V_slow = V(1); %no differentiation of slow vs. fast to start
%joint_ent=1.0; %unused at the moment.

    %NOEMO params 8 x 1 <numeric>
    %    ,1:  lambda           #weight for previous trial RT (autocorrelation of RT_t with RT_t-1)
    %    ,2:  epsilon          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %                              about fast vs. slow responses
    %    ,3:  alphaG           #learning rate for positive prediction errors (approach)
    %    ,4:  alphaN           #learning rate for negative prediction errors (avoid)
    %    ,5:  K                #baseline response speed (person mean RT?)
    %    ,6:  scale            #nu: going for the gold (modulating RT toward highest payoff)
    %    ,7:  rho         #rho parameter: weight for expected reward of fast versus slow

    lambda = params(1);
    epsilon = params(2);
    alphaG =  params(3);
    alphaN = params(4);
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(5);
    scale = params(6);
    sticky_decay = -1; %not relevant
    decay = 1;  % decay counts for beta distribution 1= nodecay
    rho = params(7);
    

    %STICKY MODEL -- OMITTING FOR NOW
    %NOEMOSTICKY params 8 x 1 <numeric>
    %    ,1:  lambda           #weight for sticky choice (weight influencing prior RTs' effect on
    %                              current RT. This is a decaying function of RT history; see sticky_decay)
    %    ,2:  epsilon          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %                              about fast vs. slow responses
    %    ,3:  alphaG           #learning rate for positive prediction errors (approach)
    %    ,4:  alphaN           #learning rate for negative prediction errors (avoid)
    %    ,5:  K                #baseline response speed (person mean RT?)
    %    ,6:  sticky_decay          #d: decay parameter influencing the degree to which prior RTs continue to affect current RTs
    %    ,7:  rho         #rho parameter: weight for expected reward of fast versus slow

%     lambda = params(1);
%     epsilon = params(2);
%     alphaG =  params(3);
%     alphaN = params(4);
%     alphaV =  0.1; % just set this to avoid degeneracy
%     K = params(5);
%     scale = -1; %going for the gold not relevant for sticky choice
%     sticky_decay = params(6);
%     decay = 1;  % decay counts for beta distribution 1= nodecay
%     rho = params(7);

% fprintf('params: k=%.2f; lambda=%.2f; nu=%.2f; alphaG=%.2f; alphaN=%.2f; rho=%.2f; epsilon=%.2f\n', ...
%     K, lambda, scale, alphaG, alphaN, rho, epsilon);

Noise=0;

%initialize algorithm parameters
exp=0;
explore=0;
mean_short = 0.5;
mean_long = 0.5;
RT_locavg = RTpred(1); % set local/learned avg on first trial..

alph_long=1.01; beta_long=1.01; % init counters and beta distribution hyperparams..
alpha_short=1.01; beta_short=1.01;

sticky      = 0;        % initialize sticky choice

if scale == -1
    %sticky choice model
    %not sure why these are zeroed at initialization for sticky choice
    RT_last2    = 0; % just for init
    RT_last3    = 0; % just for init
else
    %regular TC model including nu (going for gold)
    RT_last2    = RTpred(1); % just for init
    RT_last3    = RTpred(1); % just for init
end
bestRT      = RTpred(1);   % just for init

var_short = alpha_short*beta_short/(((alpha_short+beta_short)^2)*(alpha_short+beta_short+1));
var_long = alph_long*beta_long/(((alph_long+beta_long)^2)*(alph_long+beta_long+1));

%initialize return values for PEs, mean for short and long, etc.
ret.rtpred = NaN(ntrials, 1);
ret.rpe = NaN(ntrials, 1);
ret.explore = NaN(ntrials, 1);
ret.sdShort = NaN(ntrials, 1);
ret.sdLong = NaN(ntrials, 1);
ret.meanShort = NaN(ntrials, 1);
ret.meanLong = NaN(ntrials, 1);
ret.go = NaN(ntrials, 1);
ret.noGo = NaN(ntrials, 1);
ret.V = V;
ret.cond = cond; %add contingency to return values

%iterate over trial 2..n
for i = 2:ntrials
    explore_last = explore;
    means_last = mean_short;
    meanl_last = mean_long;
    vars_last = var_short;
    varl_last = var_long;
    goldlast = scale*(bestRT - priors.AvgRT); %just use midpoint %scale*(bestRT-RT_avg);
    bestRT_last=bestRT;
    %avg_RT_last = avg_RT;
    
    V(i) = V(i-1) +alphaV*(Reward(i-1) - V(i-1)); % update critic value
    
    rew_max = max(Reward(1:i-1)); % max reward received in block thus far -- used for v scaling v[RT_best - RT_locavg]
    rew_std = std(Reward(1:i-1)); % stddev of rewards in block thus far
    
    if Reward(i-1) > V(i-1) && Reward(i-1) >= (rew_max-rew_std), 
        % save Rt corresponding to most recent reward within one sd of max
        bestRT=RTpred(i-1);
    end;
    
    % Carry forward Go and NoGo terms unless updated by PE
    Go(i)   = Go(i-1); 
    NoGo(i) = NoGo(i-1);
    
    %CALCULATE SPEED-UP OR SLOW-DOWN OF RT DUE TO PPE OR NPE
    if(Reward(i-1) > V(i-1))
        %if obtained reward was better than expected, speed up (scaled by alphaG)
        Go(i) = Go(i-1) + alphaG*(Reward(i-1) - V(i-1));
    elseif(Reward(i-1) <= V(i-1))
        %if obtained reward was worse than expected, slow down (scaled by alphaN)
        NoGo(i) = NoGo(i-1) + alphaN*(V(i-1) - Reward(i-1));
    end
       
    %CALCULATE STRATEGIC EXPLORATION
    if(RTpred(i-1) > RT_locavg) % last response was slow/long
            
        if(Reward(i-1)> V(i-1))
            alph_long = alph_long +1; % increment count for beta distribution
        else
            beta_long = beta_long+1;
        end
        
        alph_long=decay*alph_long; % if decay <1 then this decays counts, making beta dists less confident
        beta_long=decay*beta_long;
        alpha_short=decay*alpha_short;
        beta_short=decay*beta_short;
        
        % these are mode and variances of beta dists
        var_short = alpha_short*beta_short/(((alpha_short+beta_short)^2)*(alpha_short+beta_short+1));
        var_long = alph_long*beta_long/(((alph_long+beta_long)^2)*(alph_long+beta_long+1));
        mode_long = (alph_long -1) / (alph_long + beta_long -2);
        mode_short = (alpha_short -1) / (alpha_short + beta_short -2);
        mean_long = (alph_long) / (alph_long + beta_long);
        mean_short = (alpha_short) / (alpha_short + beta_short);
        
        explore = - epsilon*(sqrt(var_short) - sqrt(var_long));  % speed up if more uncertain about fast responses
        
    elseif (RTpred(i-1)<= RT_locavg)  % last resp was fast/short
        
        % only update rew statistics if subject actually responded
        % non-response is counted as 0 in e-prime version
        if(RTpred(i-1) > 0)
            
            if(Reward(i-1)> V(i-1))
                alpha_short = alpha_short +1;
            else
                beta_short = beta_short +1;
            end
            alph_long=decay*alph_long;
            beta_long=decay*beta_long;
            alpha_short=decay*alpha_short;
            beta_short=decay*beta_short;
            
            % mode and variances of beta distribution
            var_short = alpha_short*beta_short/(((alpha_short+beta_short)^2)*(alpha_short+beta_short+1));
            var_long = alph_long*beta_long/(((alph_long+beta_long)^2)*(alph_long+beta_long+1));
            mode_long = (alph_long -1) / (alph_long + beta_long -2);
            mode_short = (alpha_short -1) / (alpha_short + beta_short -2);
            mean_long = (alph_long) / (alph_long + beta_long);
            mean_short = (alpha_short) / (alpha_short + beta_short);
            
            explore = + epsilon*(sqrt(var_long) - sqrt(var_short));
            
        end;
    end;

    % reset if already explored in this direction last trial
    % (see supplement of Frank et al 09)
    if i > 2
        if RTpred(i-1) < RTpred(i-2) && explore < 0
            explore = 0;
        elseif RTpred(i-1) > RTpred(i-2) && explore > 0
            explore=0;
        end
    end
    
    % if last trial there was no response, use trial before that for updating RT avg and autocorrelation effects (otherwise counted as 0)
    %if RTpred(i-1)==0, RTpred(i-1) = RT_last2; end; 
    
    RT_locavg = RT_locavg + alphaV*(RTpred(i-1)-RT_locavg); %update average RT locally...
    
    if scale == -1
        %sticky model: scale effect of prior RTs (decayed) on current RT by lambda
        %model does not include going for the gold (scale) update.
        RTpred(i) = K + lambda*sticky - Go(i) + NoGo(i) + explore ...
            + rho*(mean_long-mean_short) + Noise*(rand-0.5);
    else
        RTpred(i) = K + lambda*RTpred(i-1) - Go(i) + NoGo(i) + explore ...
            + rho*(mean_long-mean_short) + scale*(bestRT - priors.AvgRT) + Noise*(rand-0.5);
    end
       
    %make sure agent doesn't sample outside of valid interval
    if RTpred(i) > rtbounds(2)
        RTpred(i) = rtbounds(2);
    elseif RTpred(i) < rtbounds(1)
        RTpred(i) = rtbounds(1);
    end

    rtdiv = round(RTpred(i)/10);       
    if rtdiv < 1
      rtdiv = 1;
    elseif rtdiv > round(rtbounds(2)/10)
      rtdiv = floor(rtbounds(2)/10);
    end							 

    %Enact predicted RT to obtain reward
    if usestruct
        [Reward(i), ev(i), cstruct] = getNextRew(rtdiv, cstruct);
    else
        [Reward(i), ev(i)] = RewFunction(RTpred(i), cond, userngseed);
    end
    
    ret.rew(i-1) = Reward(i-1);             %rewards
    ret.rpe(i-1) = Reward(i-1) - V(i-1);    %prediction error
    ret.explore(i-1) = explore_last;           %explore product: epsilon * (sd diff [slow - fast])
    ret.sdShort(i-1) = sqrt(vars_last);     %sd of short RT
    ret.sdLong(i-1) = sqrt(varl_last);      %sd of long RT
    ret.meanShort(i-1) = means_last;        %mean of short RT
    ret.meanLong(i-1) = meanl_last;         %mean of long RT
    ret.go(i-1) = Go(i-1);                  %mean of Gos
    ret.noGo(i-1) = NoGo(i-1);              %mean of NoGos
    ret.V(i-1) = V(i-1);                    %expected value
    ret.rho(i-1) = rho*(meanl_last - means_last);           %shift based on mean difference between fast and slow
    ret.gold(i-1) = goldlast;
    ret.bestRT(i-1) = bestRT_last;
    %ret.avg_RT(i-1) = avg_RT_last;
    
    if (i == ntrials)
        %If this is the last trial, add return values for last trial
        ret.rew(i) = Reward(i);             %rewards
        ret.rpe(i) = Reward(i) - V(i);     %prediction error
        ret.explore(i) = explore;                  %explore product: epsilon * (sd diff [slow - fast])
        ret.sdShort(i) = sqrt(var_short);       %sd of short RT
        ret.sdLong(i) = sqrt(var_long);         %sd of long RT
        ret.meanShort(i) = mean_short;          %mean of short RT
        ret.meanLong(i) = mean_long;            %mean of long RT
        ret.go(i) = Go(i);                     %mean of Gos
        ret.noGo(i) = NoGo(i);                 %mean of NoGos
        ret.V(i) = V(i);                       %expected value
        ret.rho(i) = rho*(mean_long - mean_short);
        ret.gold(i) = scale*(bestRT-RT_locavg);
        ret.bestRT(i) = bestRT;
        %ret.avg_RT(i) = avg_RT;
    end
    
    trial_plots = 0;
    
    if trial_plots==1
        figure(1); clf;
        set(gca,'FontSize',18);
        %subplot(3,2,1);
        title('Choice history');
        %plot(tvec,v_func);
        scatter(RTpred(1:i-1),Reward(1:i-1)); axis([1 500 0 350]);
        hold on;
        plot(RTpred(i-1),Reward(i-1),'r*','MarkerSize',20);  axis([1 500 0 350]);
        hold off;
        %subplot(3,2,2)
        %title('Learned value');
        %plot(tvec,v_func); xlim([-1 ntimesteps+1]);
        %ylabel('expected value')
        %subplot(3,2,3);
        
        drawnow update;
        
    end
    
end
cost = -sum(ev);
ret.rtpred = RTpred; %copy predicted RTs to return struct (e.g., used for plotting)
ret.ev = ev;
