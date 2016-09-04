function [RTpred, ret]=TC_Alg(RTobs, Reward, params, priors, avg_RT, rewFunc, emo, model)
% Time clock R-L algorithm developed by Michael Frank
%
%
% inputs:
%    RTobs:    vector of observed reaction times
%    Reward:   vector of obtained rewards (points)
%    params:   vector of model parameters used to fit data
%    avg_RT:   scalar of average reaction time (across trials from all blocks)
%    rewFunc:  reward contingency, 1=CEV; 2=CEVR; 3=DEV; 4=IEV
%    emo:      emotion, 1=happy; 2=fear; 3=scrambled
%
%
% RTobs is a t x 1 column vector of reaction times for t trials.
% Reward is a t x 1 column vector of rewards obtained for t trials.
% params is an 8 x 1 column vector of model parameters to be used to fit behavior.
%

global emoNames;

if nargin < 7, emo = -1; end %no emotion condition
if nargin < 8, model = 'noemo'; end %if not specified, use the model where no parameters vary by emotion.

numTrials = length(RTobs);

RTpred = NaN(numTrials, 1); %vector of predicted RTs
V = NaN(numTrials, 1); %state-value function (expected value)
Go = NaN(numTrials, 1);
NoGo = NaN(numTrials, 1);

% can't predict first choice due to learning so just set it to actual subject RT
% and then predict starting trial 2
RTpred(1) = RTobs(1)*ones(1,1);

V(1) = priors.V; %initialize expected value for first trial to prior (possibly from previous run)
Go(1) = priors.Go; %initialize Go for first trial
NoGo(1) = priors.NoGo; %initialize NoGo for first trial

%V_fast = V(1); V_slow = V(1); %no differentiation of slow vs. fast to start
%joint_ent=1.0; %unused at the moment.

if strcmp(model, 'noemo') || strcmp(model, 'noemo_scram')
    %NOEMO params 8 x 1 <numeric>
    %    ,1:  lambda           #weight for previous trial RT (autocorrelation of RT_t with RT_t-1)
    %    ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %                              about fast vs. slow responses
    %    ,3:  alpha1           #learning rate for positive prediction errors (approach)
    %    ,4:  alpha2           #learning rate for negative prediction errors (avoid)
    %    ,5:  K                #baseline response speed (person mean RT?)
    %    ,6:  scale            #nu: going for the gold (modulating RT toward highest payoff)
    %    ,7:  exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %    ,8:  meandiff         #rho parameter: weight for expected reward of fast versus slow

    lambda = params(1);
    explore = params(2);
    alpha1 =  params(3);
    alpha2 = params(4);
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(5);
    scale = params(6);
    sticky_decay = -1; %not relevant
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(7); % param for alternative exp models of rt swings
    meandiff = params(8);
    
elseif strcmp(model, 'noemosticky') || strcmp(model, 'noemosticky_scram')
    %NOEMOSTICKY params 8 x 1 <numeric>
    %    ,1:  lambda           #weight for sticky choice (weight influencing prior RTs' effect on
    %                              current RT. This is a decaying function of RT history; see sticky_decay)
    %    ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %                              about fast vs. slow responses
    %    ,3:  alpha1           #learning rate for positive prediction errors (approach)
    %    ,4:  alpha2           #learning rate for negative prediction errors (avoid)
    %    ,5:  K                #baseline response speed (person mean RT?)
    %    ,6:  sticky_decay          #d: decay parameter influencing the degree to which prior RTs continue to affect current RTs
    %    ,7:  exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %    ,8:  meandiff         #rho parameter: weight for expected reward of fast versus slow

    lambda = params(1);
    explore = params(2);
    alpha1 =  params(3);
    alpha2 = params(4);
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(5);
    scale = -1; %going for the gold not relevant for sticky choice
    sticky_decay = params(6);
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(7); % param for alternative exp models of rt swings
    meandiff = params(8);
    
elseif strcmp(model, 'emoexplore')
    %EMOEXPLORE params 10 x 1 <numeric>
    %
    %   ,1:  lambda           #weight for previous trial RT (autocorrelation of RT_t with RT_t-1)
    %   ,2:  explore_scram    #epsilon parameter for scrambled: how much should RT be modulated by greater relative uncertainty
    %   ,3:  explore_fear     #epsilon parameter for fearful: how much should RT be modulated by greater relative uncertainty
    %   ,4:  explore_happy    #epsilon parameter for happy: how much should RT be modulated by greater relative uncertainty
    %   ,5:  alpha1           #learning rate for positive prediction errors (approach)
    %   ,6:  alpha2           #learning rate for negative prediction errors (avoid)
    %   ,7:  K                #baseline response speed (person mean RT?)
    %   ,8:  scale            #nu: going for the gold (modulating RT toward highest payoff)
    %   ,9:  exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %   ,10: meandiff         #rho parameter: weight for expected reward of fast versus slow
    
    lambda = params(1);
    if strcmp(emoNames{emo}, 'happy')
        explore=params(4);
    elseif strcmp(emoNames{emo}, 'fear')
        explore=params(3);
    elseif strcmp(emoNames{emo}, 'scram')
        explore = params(2);
    end
    alpha1 =  params(5);
    alpha2 = params(6);
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(7);
    scale = params(8);
    sticky_decay = -1; %not relevant
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(9); % param for alternative exp models of rt swings
    meandiff = params(10);
    
elseif strcmp(model, 'emoexploresticky')
    %EMOEXPLORESTICKY params 10 x 1 <numeric>
    %
    %   ,1:  lambda           #weight for sticky choice (weight influencing prior RTs' effect on
    %                              current RT. This is a decaying function of RT history; see sticky_decay)
    %   ,2:  explore_scram    #epsilon parameter for scrambled: how much should RT be modulated by greater relative uncertainty
    %   ,3:  explore_fear     #epsilon parameter for fearful: how much should RT be modulated by greater relative uncertainty
    %   ,4:  explore_happy    #epsilon parameter for happy: how much should RT be modulated by greater relative uncertainty
    %   ,5:  alpha1           #learning rate for positive prediction errors (approach)
    %   ,6:  alpha2           #learning rate for negative prediction errors (avoid)
    %   ,7:  K                #baseline response speed (person mean RT?)
    %   ,8:  sticky_decay     #d: decay parameter influencing the degree to which prior RTs continue to affect current RTs
    %   ,9:  exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %   ,10: meandiff         #rho parameter: weight for expected reward of fast versus slow
    
    lambda = params(1);
    if strcmp(emoNames{emo}, 'happy')
        explore=params(4);
    elseif strcmp(emoNames{emo}, 'fear')
        explore=params(3);
    elseif strcmp(emoNames{emo}, 'scram')
        explore = params(2);
    end
    alpha1 =  params(5);
    alpha2 = params(6);
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(7);
    scale = -1; %going for the gold not relevant for sticky choice
    sticky_decay = params(8);
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(9); % param for alternative exp models of rt swings
    meandiff = params(10);
    
elseif strcmp(model, 'emogonogo')
    %EMOGONOGO: model parameter initialization
    %Params 12 x 1 <numeric>
    %   ,1:  lambda           #weight for previous trial RT (autocorrelation of RT_t with RT_t-1)
    %   ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %   ,3:  alpha1_scram     #learning rate for scrambled positive prediction errors (approach)
    %   ,4:  alpha1_fear      #learning rate for fear positive prediction errors (approach)
    %   ,5:  alpha1_happy     #learning rate for happy positive prediction errors (approach)
    %   ,6:  alpha2_scram     #learning rate for scrambled negative prediction errors (avoid)
    %   ,7:  alpha2_fear      #learning rate for fear negative prediction errors (avoid)
    %   ,8:  alpha2_happy     #learning rate for happy negative prediction errors (avoid)
    %   ,9:  K                #baseline response speed (person mean RT?)
    %   ,10: scale            #nu: going for the gold (modulating RT toward highest payoff)
    %   ,11: exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %   ,12: meandiff         #rho parameter: weight for expected reward of fast versus slow
    
    lambda = params(1);
    explore = params(2);
    if strcmp(emoNames{emo}, 'happy')
        alpha1=params(5);
        alpha2=params(8);
    elseif strcmp(emoNames{emo}, 'fear')
        alpha1=params(4);
        alpha2=params(7);
    elseif strcmp(emoNames{emo}, 'scram')
        alpha1=params(3);
        alpha2=params(6);
    end
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(9);
    scale = params(10);
    sticky_decay = -1; %not relevant
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(11); % param for alternative exp models of rt swings
    meandiff = params(12);
    
    
elseif strcmp(model, 'emogonogosticky')
    %EMOGONOGOSTICKY: model parameter initialization
    %Params 12 x 1 <numeric>
    %   ,1:  lambda           #weight for sticky choice (weight influencing prior RTs' effect on
    %                              current RT. This is a decaying function of RT history; see sticky_decay)
    %   ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %   ,3:  alpha1_scram     #learning rate for scrambled positive prediction errors (approach)
    %   ,4:  alpha1_fear      #learning rate for fear positive prediction errors (approach)
    %   ,5:  alpha1_happy     #learning rate for happy positive prediction errors (approach)
    %   ,6:  alpha2_scram     #learning rate for scrambled negative prediction errors (avoid)
    %   ,7:  alpha2_fear      #learning rate for fear negative prediction errors (avoid)
    %   ,8:  alpha2_happy     #learning rate for happy negative prediction errors (avoid)
    %   ,9:  K                #baseline response speed (person mean RT?)
    %   ,10: sticky_decay     #d: decay parameter influencing the degree to which prior RTs continue to affect current RTs
    %   ,11: exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %   ,12: meandiff         #rho parameter: weight for expected reward of fast versus slow
    
    lambda = params(1);
    explore = params(2);
    if strcmp(emoNames{emo}, 'happy')
        alpha1=params(5);
        alpha2=params(8);
    elseif strcmp(emoNames{emo}, 'fear')
        alpha1=params(4);
        alpha2=params(7);
    elseif strcmp(emoNames{emo}, 'scram')
        alpha1=params(3);
        alpha2=params(6);
    end
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(9);
    scale = -1; %not relevant
    sticky_decay = params(10);
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(11); % param for alternative exp models of rt swings
    meandiff = params(12);
    
elseif strcmp(model, 'emonogo')
    %EMONOGO: initial fit of emogonogo suggested that the go parameters (PPE) had little variation.
    %         Thus, allow just no go learning rate to vary by condition (scrambled, fear, happy)
    %Params 10 x 1 <numeric>
    %   ,1:  lambda           #weight for previous trial RT (autocorrelation of RT_t with RT_t-1)
    %   ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %   ,3:  alpha1           #learning rate for positive prediction errors (approach)
    %   ,4:  alpha2_scram     #learning rate for scrambled negative prediction errors (avoid)
    %   ,5:  alpha2_fear      #learning rate for fear negative prediction errors (avoid)
    %   ,6:  alpha2_happy     #learning rate for happy negative prediction errors (avoid)
    %   ,7:  K                #baseline response speed (person mean RT?)
    %   ,8:  scale            #nu: going for the gold (modulating RT toward highest payoff)
    %   ,9:  exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %   ,10: meandiff         #rho parameter: weight for expected reward of fast versus slow
    
    lambda = params(1);
    explore = params(2);
    alpha1 = params(3);
    if strcmp(emoNames{emo}, 'happy')
        alpha2=params(6);
    elseif strcmp(emoNames{emo}, 'fear')
        alpha2=params(5);
    elseif strcmp(emoNames{emo}, 'scram')
        alpha2=params(4);
    end
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(7);
    scale = params(8);
    sticky_decay = -1; %not relevant
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(9); % param for alternative exp models of rt swings
    meandiff = params(10);
    
elseif strcmp(model, 'emonogosticky')
    %EMONOGOSTICKY: initial fit of emogonogo suggested that the go parameters (PPE) had little variation.
    %         Thus, allow just no go learning rate to vary by condition (scrambled, fear, happy)
    %Params 10 x 1 <numeric>
    %   ,1:  lambda           #weight for sticky choice (weight influencing prior RTs' effect on
    %                              current RT. This is a decaying function of RT history; see sticky_decay)
    %   ,2:  explore          #epsilon parameter: how much should RT be modulated by greater relative uncertainty
    %   ,3:  alpha1           #learning rate for positive prediction errors (approach)
    %   ,4:  alpha2_scram     #learning rate for scrambled negative prediction errors (avoid)
    %   ,5:  alpha2_fear      #learning rate for fear negative prediction errors (avoid)
    %   ,6:  alpha2_happy     #learning rate for happy negative prediction errors (avoid)
    %   ,7:  K                #baseline response speed (person mean RT?)
    %   ,8:  sticky_decay          #d: decay parameter influencing the degree to which prior RTs continue to affect current RTs
    %   ,9:  exp_alt          #alternative exponential models for RT swings (not sure of its use yet)
    %   ,10: meandiff         #rho parameter: weight for expected reward of fast versus slow
    
    lambda = params(1);
    explore = params(2);
    alpha1 = params(3);
    if strcmp(emoNames{emo}, 'happy')
        alpha2=params(6);
    elseif strcmp(emoNames{emo}, 'fear')
        alpha2=params(5);
    elseif strcmp(emoNames{emo}, 'scram')
        alpha2=params(4);
    end
    alphaV =  0.1; % just set this to avoid degeneracy
    K = params(7);
    scale = -1; %not relevant
    sticky_decay = params(8);
    decay = 1;  % decay counts for beta distribution 1= nodecay
    exp_alt = params(9); % param for alternative exp models of rt swings
    meandiff = params(10);
end

%vary params 

Q = 0;

Noise=0;

dist_type = 'beta';

mean_s = 0; mean_f = 0;
%initial variances of fast/slow resps for kalman filter
%just init to rewvar so initial lr = 0.5
vars = (std(Reward))^2;varf = (std(Reward))^2;

% GENERATIVE model just pick some params to generate data
global generative;
if (generative ==1)
    lambda = 0.2;
    explore = 3000;
    alpha1 = .3;
    alpha2 = .3;
    K = 1500;
    exp_alt =500;
    scale = .25;
    meandiff = 1000;
    if strcmp(dist_type,'Gauss')
        meandiff=20;
        explore = 10;
    end
    Noise=2000;
end

%initialize algorithm parameters
exp=0;
exp1=0;
exp1a=0;
lose_switch = 0;
regress = 0;
mean_short = 0.5;
mean_long = 0.5;
RT_avg = avg_RT; % set avg on first trial..

alph_long=1.01; b_long=1.01; % init counters and beta distribution hyperparams..
alph_short=1.01; b_short=1.01;
cnt_short=0; cnt_long=0;
cnt_speed=0; cnt_slow=0;

RT_new      = RTobs(1); % just for init
RT_last     = RTobs(1); % just for init
sticky      = 0;        % initialize sticky choice

if scale == -1
    %sticky choice model
    %not sure why these are zeroed at initialization for sticky choice
    RT_last2    = 0; % just for init
    RT_last3    = 0; % just for init
else
    %regular TC model including nu (going for gold)
    RT_last2    = RTobs(1); % just for init
    RT_last3    = RTobs(1); % just for init
end
bestRT      = avg_RT;   % just for init

var_short = alph_short*b_short/(((alph_short+b_short)^2)*(alph_short+b_short+1));
var_long = alph_long*b_long/(((alph_long+b_long)^2)*(alph_long+b_long+1));

%initialize return values for PEs, mean for short and long, etc.
ret.rtobs = RTobs;
ret.rtpred = NaN(numTrials, 1);
ret.rpe = NaN(numTrials, 1);
ret.explore = NaN(numTrials, 1);
ret.sdShort = NaN(numTrials, 1);
ret.sdLong = NaN(numTrials, 1);
ret.meanShort = NaN(numTrials, 1);
ret.meanLong = NaN(numTrials, 1);
ret.go = NaN(numTrials, 1);
ret.noGo = NaN(numTrials, 1);
ret.ev = V;
ret.rewFunc = rewFunc; %add contingency to return values
ret.emo = emo;         %add emotion to return values

%iterate over trial 2..n
for trial = 2:numTrials
    lasttrial = trial-1;
    
    exp1_last = exp1; exp_last = exp; exp1a_last = exp1a;
    means_last = mean_short;
    meanl_last = mean_long;
    vars_last = var_short;
    varl_last = var_long;
    meandiff_last = meandiff;
    goldlast = scale*(bestRT-avg_RT);
    bestRT_last=bestRT;
    avg_RT_last = avg_RT;
    
    if strcmp(dist_type,'Gauss')
        % add process noise to kalman variances (only for kalman filter model)
        vars = vars+ Q; varf = varf + Q;
    end
    
    if(generative==1) % if generating responses make last rt the prev predicted rt
        if trial > 2, RT_last2= RT_last; end
        RT_last = RT_new;
    else
        RT_last = RTobs(lasttrial);
        if scale == -1, sticky = RT_last + sticky_decay*sticky; end %update sticky choice if used
        if trial > 2, RT_last2 = RTobs(trial-2); end
        if trial > 3, RT_last3 = RTobs(trial-3); end
    end
    
    mom = (RT_last - RT_last2); % momentum
    if mom==0, mom=1; end
    
    Rew_last = Reward(lasttrial);
    if (generative==1)
        Rew_last = RewFunction(RT_last, rewFunc); % calculate reward if model generating own rt's
    end
    
    if RT_last > RT_last2   % Reverse-momentum model
        cnt_slow = cnt_slow + 1; cnt_speed = 0; % count number of responses slower than previous
    else
        cnt_speed = cnt_speed + 1; cnt_slow = 0;
    end
    
    V_last = V(lasttrial);
    V_new = V_last +alphaV*(Rew_last - V_last); % update critic value
    
    rew_max = max(Reward(1:lasttrial)); % max reward received in block thus far -- used for v scaling v[RT_best - RT_avg]
    rew_std = std(Reward(1:lasttrial)); % stddev of rewards in block thus far
    
    if Rew_last > V_last && Rew_last>= (rew_max-rew_std), 
        % save Rt corresponding to most recent reward within one sd of max
        bestRT=RT_last;
    end;
    
    Go_last = Go(lasttrial);
    NoGo_last = NoGo(lasttrial);
    Go_new = Go_last; NoGo_new = NoGo_last; %  unless updated below by PE
    
    %process speed-up or slow-down of RT due to PPE or NPE
    if(Rew_last> V_last)
        %if obtained reward was better than expected, speed up (scaled by alpha1)
        Go_new = Go_last + alpha1*(Rew_last - V_last);
    elseif(Rew_last <= V_last)
        %if obtained reward was worse than expected, slow down (scaled by alpha1)
        NoGo_new = NoGo_last + alpha2*(V_last - Rew_last);
    end
    
    if(RT_last > RT_avg) % last response was slow/long
        cnt_long = cnt_long +1; cnt_short =0; % for sutton exp
        % bonus control model
        
        regress = -exp_alt;  % for simple oscillation regression to
        % mean explore control model
        
        if strcmp(dist_type,'beta')
            
            if(Rew_last> V_last)
                alph_long = alph_long +1; % increment count for beta distribution
            else
                b_long = b_long+1;
                lose_switch = -exp_alt; % if was slow go fast after neg PE
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
            
        elseif strcmp(dist_type,'Gauss')
            
            rewvar = (std(Reward))^2;
            
            alphaKs = vars/(vars+rewvar); % Kalman gain for slow responses
            vars = (1 - alphaKs)*vars; % Kalaman variance for slow resps
            
            mean_s = mean_s + alphaKs*((Rew_last - 0*V_last) - mean_s); ...
                % kalman mean
            
            mean_long = mean_s; mean_short=mean_f;
            var_short =varf; var_long = vars;
            
            exp1 = - explore*(sqrt(varf) - sqrt(vars));  % using kalman filter gaussian distributions.
            
        end
        
        if RT_last < RT_last2 && exp1 < 0, exp1=0; %-exp1; % reset if
            %already explored in this direction last trial (see supplement of  Frank et al 09)
        elseif RT_last>RT_last2 && exp1 > 0, exp1=0;
        end;
        
    elseif (RT_last<= RT_avg)  % last resp was fast/short
        
         % only update rew statistics if subject actually responded
         % non-response is counted as 0 in e-prime version
        if(RT_last > 0)
            
            cnt_short = cnt_short +1; cnt_long =0; % for sutton exp
            % bonus control model
            
            regress = +exp_alt; %  regress to mean control model
            
            if strcmp(dist_type,'beta')
                
                if(Rew_last> V_last)
                    alph_short = alph_short +1;
                else
                    b_short = b_short +1;
                    lose_switch = exp_alt; % if was fast, slow down after neg PE (lose switch control model)
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
                
            elseif strcmp(dist_type,'Gauss') % for kalman filter model
                
                rewvar = (std(Reward))^2;
                alphaKf = varf / (varf +rewvar);
                varf = (1 - alphaKf)*varf;
                mean_f = mean_f + alphaKf*((Rew_last - 0*V_last) - mean_f);
                mean_short = mean_f; mean_long = mean_s;
                var_short =varf; var_long = vars;
                
                exp1 = + explore*(sqrt(vars) - sqrt(varf));  % using kalman filter normal distributions.
            end
            
            % reset if already explored in this direction last trial
            % (see supplement of Frank et al 09)
            if RT_last < RT_last2 && exp1 < 0
                exp1 = 0;
            elseif RT_last > RT_last2 && exp1 > 0
                exp1=0;
            end;
            
        end;
    end;
    
    revmom=0; % for reverse momentum control model
    if cnt_speed > scale, revmom=exp_alt*cnt_speed;
    elseif cnt_slow > scale, revmom = -exp_alt*cnt_slow;
    end
    
    exp = exp_alt*(sqrt(cnt_short)-sqrt(cnt_long)); % sutton
    % exploration bonus model, for control model in supplement
    
    if RT_last==0, RT_last = RT_last2; end; %% if last trial there
    %% was no response, use trial before that for updating RT avg and autocorrelation effects (otherwise counted as 0)
    
    RT_avg = RT_avg + alphaV*(RT_last-RT_avg); %update average RT locally...
    
    if scale == -1
        %sticky model: scale effect of prior RTs (decayed) on current RT by lambda
        %model does not include going for the gold (scale) update.
        RT_new = K + lambda*sticky - Go_new + NoGo_new  +exp1 + 0*regress ...
            + meandiff*(mean_long-mean_short) + Noise*(rand-0.5);
    else
        RT_new = K + lambda*RT_last - Go_new + NoGo_new  +exp1 + 0*regress ...
            + meandiff*(mean_long-mean_short) + scale*(bestRT-avg_RT)+ Noise*(rand-0.5);
    end
    
    if RTobs(trial)==0, RT_new = 0; end; %% don't try to predict response failures, which are counted as 0 in e-prime
    
    rtvar = (std(RTobs))^2;
    
    if strcmp(dist_type,'Gauss')
        alphaK = varK/(varK+rtvar); % alphaK = kalman gain;
        varK = (1 - alphaK)*varK;
    end
    
    RTpred(trial) = RT_new; % add new RT pred to vector
    V(trial) = V_new;
    Go(trial) = Go_new;
    NoGo(trial) = NoGo_new;
    
    ret.rpe(trial-1) = Rew_last - V_last;       %prediction error
    ret.explore(trial-1) = exp1_last;           %explore product: epsilon * (sd diff [slow - fast])
    ret.sdShort(trial-1) = sqrt(vars_last);     %sd of short RT
    ret.sdLong(trial-1) = sqrt(varl_last);      %sd of long RT
    ret.meanShort(trial-1) = means_last;        %mean of short RT
    ret.meanLong(trial-1) = meanl_last;         %mean of long RT
    ret.go(trial-1) = Go_last;                  %mean of Gos
    ret.noGo(trial-1) = NoGo_last;              %mean of NoGos
    ret.ev(trial-1) = V_last;                   %expected value
    ret.rho(trial-1) = meandiff_last*(meanl_last - means_last);           %shift based on mean difference between fast and slow
    ret.gold(trial-1) = goldlast;
    ret.bestRT(trial-1) = bestRT_last;
    ret.avg_RT(trial-1) = avg_RT_last;
    
    if (trial == numTrials)
        %If this is the last trial, add return values for last trial
        ret.rpe(trial) = Reward(trial) - V_new;     %prediction error
        ret.explore(trial) = exp1;                  %explore product: epsilon * (sd diff [slow - fast])
        ret.sdShort(trial) = sqrt(var_short);       %sd of short RT
        ret.sdLong(trial) = sqrt(var_long);         %sd of long RT
        ret.meanShort(trial) = mean_short;          %mean of short RT
        ret.meanLong(trial) = mean_long;            %mean of long RT
        ret.go(trial) = Go_new;                     %mean of Gos
        ret.noGo(trial) = NoGo_new;                 %mean of NoGos
        ret.ev(trial) = V_new;                      %expected value
        ret.rho(trial) = meandiff*(mean_long - mean_short);
        ret.gold(trial) = scale*(bestRT-avg_RT);
        ret.bestRT(trial) = bestRT;
        ret.avg_RT(trial) = avg_RT;
    end
    
end

ret.rtpred = RTpred; %copy predicted RTs to return struct (e.g., used for plotting)
