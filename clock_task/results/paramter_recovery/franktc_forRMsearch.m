function [cost, RTpred, ret]=franktc_forRMsearch(params, priors,rt_obs, rew_obs,cond, rngseeds, ntrials, rtbounds)
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
%

% if nargin < 4, rngseeds=31; end
% if nargin < 5, ntrials=50; end
% if nargin < 6, rtbounds=[0 5000]; end

%fprintf('Fitting cond %s\n', cond);

global rew_rng_state;

%initialize states for two repeatable random number generators using different seeds
rew_rng_seed=rngseeds(1);

rng(rew_rng_seed);
rew_rng_state=rng;

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

rt_obs(1) = 2500; %start in middle
RTpred(1) = rt_obs(1);

% if usestruct
%     [Reward(1), ev(1), cstruct] = getNextRew(rt_obs(1)/10, cstruct);
% else
%     [Reward(1), ev(1)] = RewFunction(rt_obs(1), cond);
% end

Reward(1) = rew_obs(1); %Initialize the first reward
V(1) = priors.V; %initialize expected value for first trial to prior (possibly from previous run)
Go(1) = priors.Go; %initialize Go for first trial
NoGo(1) = priors.NoGo; %initialize NoGo for first trial

%V_fast = V(1); V_slow = V(1); %no differentiation of slow vs. fast to start
%joint_ent=1.0; %unused at the moment.

    %NOEMO params 8 x 1 <numeric>
    %    ,1:  lambda           #weight for previous trial RT (autocorrelation of RT_t with RT_t-1)
    %    ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %                              about fast vs. slow responses
    %    ,3:  alphaG           #learning rate for positive prediction errors (approach)
    %    ,4:  alphaN           #learning rate for negative prediction errors (avoid)
    %    ,5:  K                #baseline response speed (person mean RT?)
    %    ,6:  scale            #nu: going for the gold (modulating RT toward highest payoff)
    %    ,7:  meandiff         #rho parameter: weight for expected reward of fast versus slow

    lambda = params(1);
    explore = params(2);
    alphaG =  params(3);
    alphaN = params(4);
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(5);
    scale = params(6);
    sticky_decay = -1; %not relevant
    decay = 1;  % decay counts for beta distribution 1= nodecay
    meandiff = params(7);
    

    %STICKY MODEL -- OMITTING FOR NOW
    %NOEMOSTICKY params 8 x 1 <numeric>
    %    ,1:  lambda           #weight for sticky choice (weight influencing prior RTs' effect on
    %                              current RT. This is a decaying function of RT history; see sticky_decay)
    %    ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %                              about fast vs. slow responses
    %    ,3:  alphaG           #learning rate for positive prediction errors (approach)
    %    ,4:  alphaN           #learning rate for negative prediction errors (avoid)
    %    ,5:  K                #baseline response speed (person mean RT?)
    %    ,6:  sticky_decay          #d: decay parameter influencing the degree to which prior RTs continue to affect current RTs
    %    ,7:  meandiff         #rho parameter: weight for expected reward of fast versus slow

%     lambda = params(1);
%     explore = params(2);
%     alphaG =  params(3);
%     alphaN = params(4);
%     alphaV =  0.1; % just set this to avoid degeneracy
%     K = params(5);
%     scale = -1; %going for the gold not relevant for sticky choice
%     sticky_decay = params(6);
%     decay = 1;  % decay counts for beta distribution 1= nodecay
%     meandiff = params(7);

%vary params 

Noise=0;

dist_type = 'beta';

%initialize algorithm parameters
exp=0;
exp1=0;
mean_short = 0.5;
mean_long = 0.5;
RT_locavg = rt_obs(1); % set local/learned avg on first trial..

alph_long=1.01; b_long=1.01; % init counters and beta distribution hyperparams..
alph_short=1.01; b_short=1.01;

sticky      = 0;        % initialize sticky choice

if scale == -1
    %sticky choice model
    %not sure why these are zeroed at initialization for sticky choice
    RT_last2    = 0; % just for init
    RT_last3    = 0; % just for init
else
    %regular TC model including nu (going for gold)
    RT_last2    = rt_obs(1); % just for init
    RT_last3    = rt_obs(1); % just for init
end
bestRT      = rt_obs(1);   % just for init

var_short = alph_short*b_short/(((alph_short+b_short)^2)*(alph_short+b_short+1));
var_long = alph_long*b_long/(((alph_long+b_long)^2)*(alph_long+b_long+1));

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
    exp1_last = exp1;
    means_last = mean_short;
    meanl_last = mean_long;
    vars_last = var_short;
    varl_last = var_long;
    goldlast = scale*(bestRT-2500); %just use midpoint %scale*(bestRT-RT_locavg);
    bestRT_last=bestRT;
    %avg_RT_last = avg_RT;
    
    V(i) = V(i-1) +alphaV*(Reward(i-1) - V(i-1)); % update critic value
    
    rew_max = max(Reward(1:i-1)); % max reward received in block thus far -- used for v scaling v[RT_best - RT_locavg]
    rew_std = std(Reward(1:i-1)); % stddev of rewards in block thus far
    
    if Reward(i-1) > V(i-1) && Reward(i-1) >= (rew_max-rew_std), 
        % save Rt corresponding to most recent reward within one sd of max
        bestRT=rt_obs(i-1);
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
    if(rt_obs(i-1) > RT_locavg) % last response was slow/long
            
        if(Reward(i-1)> V(i-1))
            alph_long = alph_long +1; % increment count for beta distribution
        else
            b_long = b_long+1;
        end
        
        alph_long=decay*alph_long; % if decay <1 then this decays counts, making beta dists less confident
        b_long=decay*b_long;
        alph_short=decay*alph_short;
        b_short=decay*b_short;
        
        % these are mode and variances of beta dists
        var_short = alph_short*b_short/(((alph_short+b_short)^2)*(alph_short+b_short+1));
        var_long = alph_long*b_long/(((alph_long+b_long)^2)*(alph_long+b_long+1));
        mode_long = (alph_long -1) / (alph_long + b_long -2);
        mode_short = (alph_short -1) / (alph_short + b_short -2);
        mean_long = (alph_long) / (alph_long + b_long);
        mean_short = (alph_short) / (alph_short + b_short);
        
        exp1 = - explore*(sqrt(var_short) - sqrt(var_long));  % speed up if more uncertain about fast responses
        
    elseif (rt_obs(i-1)<= RT_locavg)  % last resp was fast/short
        
        % only update rew statistics if subject actually responded
        % non-response is counted as 0 in e-prime version
        if(rt_obs(i-1) > 0)
            
            if(Reward(i-1)> V(i-1))
                alph_short = alph_short +1;
            else
                b_short = b_short +1;
            end
            alph_long=decay*alph_long;
            b_long=decay*b_long;
            alph_short=decay*alph_short;
            b_short=decay*b_short;
            
            % mode and variances of beta distribution
            var_short = alph_short*b_short/(((alph_short+b_short)^2)*(alph_short+b_short+1));
            var_long = alph_long*b_long/(((alph_long+b_long)^2)*(alph_long+b_long+1));
            mode_long = (alph_long -1) / (alph_long + b_long -2);
            mode_short = (alph_short -1) / (alph_short + b_short -2);
            mean_long = (alph_long) / (alph_long + b_long);
            mean_short = (alph_short) / (alph_short + b_short);
            
            exp1 = + explore*(sqrt(var_long) - sqrt(var_short));
            
        end;
    end;

    % reset if already explored in this direction last trial
    % (see supplement of Frank et al 09)
    if i > 2
        if rt_obs(i-1) < rt_obs(i-2) && exp1 < 0
            exp1 = 0;
        elseif rt_obs(i-1) > rt_obs(i-2) && exp1 > 0
            exp1=0;
        end
    end
    
    % if last trial there was no response, use trial before that for updating RT avg and autocorrelation effects (otherwise counted as 0)
    %if RTpred(i-1)==0, RTpred(i-1) = RT_last2; end; 
    
    RT_locavg = RT_locavg + alphaV*(rt_obs(i-1)-RT_locavg); %update average RT locally...
   
    %% AD: this is the only place were you keep the RTpred
    if scale == -1
        %sticky model: scale effect of prior RTs (decayed) on current RT by lambda
        %model does not include going for the gold (scale) update.
        RTpred(i) = K + lambda*sticky - Go(i) + NoGo(i)  +exp1 ...
            + meandiff*(mean_long-mean_short) + Noise*(rand-0.5);
    else
        RTpred(i) = K + lambda*RTpred(i-1) - Go(i) + NoGo(i)  +exp1 ...
            + meandiff*(mean_long-mean_short) + scale*(bestRT-2500) + Noise*(rand-0.5); %scale*(bestRT-RT_locavg)
    end
       
    %make sure agent doesn't sample outside of valid interval
    if RTpred(i) > rtbounds(2)
        RTpred(i) = rtbounds(2);
    elseif RTpred(i) < rtbounds(1)
        RTpred(i) = rtbounds(1);
    end
   %%    
    %Enact predicted RT to obtain reward
%     if usestruct
%         %[Reward(i), ev(i), cstruct] = getNextRew(RTpred(i)/10, cstruct);
%         [Reward(i), ev(i), cstruct] = getNextRew(ceil(rt_obs(i)/10), cstruct); %Throwing an error for the lookup in getNextReward
%     else
%         [Reward(i), ev(i)] = RewFunction(rt_obs(i), cond);
%     end
    
    %Enact RT obs to obtain reward
    Reward(i) = rew_obs(i);
    
    
    ret.rew(i-1) = Reward(i-1);             %rewards
    ret.rpe(i-1) = Reward(i-1) - V(i-1);    %prediction error
    ret.explore(i-1) = exp1_last;           %explore product: epsilon * (sd diff [slow - fast])
    ret.sdShort(i-1) = sqrt(vars_last);     %sd of short RT
    ret.sdLong(i-1) = sqrt(varl_last);      %sd of long RT
    ret.meanShort(i-1) = means_last;        %mean of short RT
    ret.meanLong(i-1) = meanl_last;         %mean of long RT
    ret.go(i-1) = Go(i-1);                  %mean of Gos
    ret.noGo(i-1) = NoGo(i-1);              %mean of NoGos
    ret.V(i-1) = V(i-1);                    %expected value
    ret.rho(i-1) = meandiff*(meanl_last - means_last);           %shift based on mean difference between fast and slow
    ret.gold(i-1) = goldlast;
    ret.bestRT(i-1) = bestRT_last;
    %ret.avg_RT(i-1) = avg_RT_last;
    
    if (i == ntrials)
        %If this is the last trial, add return values for last trial
        ret.rew(i) = Reward(i);             %rewards
        ret.rpe(i) = Reward(i) - V(i);     %prediction error
        ret.explore(i) = exp1;                  %explore product: epsilon * (sd diff [slow - fast])
        ret.sdShort(i) = sqrt(var_short);       %sd of short RT
        ret.sdLong(i) = sqrt(var_long);         %sd of long RT
        ret.meanShort(i) = mean_short;          %mean of short RT
        ret.meanLong(i) = mean_long;            %mean of long RT
        ret.go(i) = Go(i);                     %mean of Gos
        ret.noGo(i) = NoGo(i);                 %mean of NoGos
        ret.V(i) = V(i);                       %expected value
        ret.rho(i) = meandiff*(mean_long - mean_short);
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
        scatter(rt_obs(1:i-1),Reward(1:i-1)); axis([1 500 0 350]);
        hold on;
        plot(rt_obs(i-1),Reward(i-1),'r*','MarkerSize',20);  axis([1 500 0 350]);
        hold off;
        %subplot(3,2,2)
        %title('Learned value');
        %plot(tvec,v_func); xlim([-1 ntimesteps+1]);
        %ylabel('expected value')
        %subplot(3,2,3);
        
        drawnow update;
        
    end
    
end
%cost = -sum(ev);
cost = sum((rt_obs-RTpred').^2)*.01; %put it on the same scale as kalman models since kalman rt_obs = 0-500 scale
ret.rtpred = RTpred; %copy predicted RTs to return struct (e.g., used for plotting)
ret.ev = ev;