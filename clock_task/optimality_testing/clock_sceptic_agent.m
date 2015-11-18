function [cost,v_it,rts,ret] = clock_sceptic_agent(params, agent, rngseeds, cond, ntrials, nbasis, ntimesteps, reversal)
% This is the primary script that attempts to solve clock contingencies using variants of the sceptic agent.
% Variants are described using the parameter agent, which determines the behavior of the forward model.
% This script is used for testing the optimality of variants of sceptic for solving clock contingencies.
%
% cond is the character string or structure of the reward contingency
% ntrials is the number of trials to run
% nbasis is the number of radial basis functions used to estimate value and uncertainty
% ntimesteps is the number of time bins used for obtaining estimates of time functions for plotting etc.
%
% Agents:
%   1) fixedLR_softmax:         fixed learning rate (alpha) for PE+ and PE-; softmax choice rule
%   2) fixedLR_egreedy:         fixed learning rate (alpha) for PE+ and PE-; epsilon greedy choice rule with rand unif exploration
%   3) fixedLR_egreedy_grw:     fixed learning rate (alpha) for PE+ and PE-; epsilon greedy choice rule with GRW exploration around rt(i)
%   4) asymfixedLR_softmax:     separate fixed learning rates (alpha and beta) for positive and negative PEs, respectively; softmax choice
%   5) kalman_softmax:          kalman learning rule (no free parameter); softmax choice over value curve
%   6) kalman_processnoise:     kalman learning rule (no free parameter); PEs enhance gain through process noise Q according to parameter omega
%   7) kalman_sigmavolatility:  kalman learning rule (no free parameter); PEs inflate posterior variance (sigma) according to phi and gamma
%   8) kalman_uv_logistic:      old kalman with explore/exploit hardmax selection according to logistic function
%   9) kalman_uv_sum:           kalman learning rule and uncertainty update; V and U are mixed by tau; softmax choice over U+V
%  10) fixedLR_kl_softmax:      fixed learning rate (alpha) for PE+ and PE-; KL discounting of value curve by PEs; softmax choice
%  11) kalman_kl_softmax:       kalman learning rule; KL discounting of value curve by PEs
%  12) kalman_processnoise_kl:  kalman learning rule; PEs enhance gain through process noise Q; KL discounting of value curve by PEs; softmax
%  13) kalman_uv_sum_kl;        kalman learning rule; V and U are mixed by tau; KL discounting of value curve by PEs; softmax choice over U+V

if nargin < 2, agent = 'fixedLR_softmax'; end
if nargin < 3, rngseeds=[98 83 66 10]; end
if nargin < 4, cond = 'DEV'; end
if nargin < 5, ntrials=100; end
if nargin < 6, nbasis = 24; end
if nargin < 7, ntimesteps=500; end
if nargin < 8, reversal = 0; end %No reversal by default

omega=0;                    %if not a PE-enhanced process noise model, remove this influence
gamma=0;                    %zero gamma and phi for all models except kalman_sigmavolatility
phi=0;
tradeoff=0;                 %parameter only applies to kalman_uv_logistic
prop_spread = params(1);    %proportion of discrete interval over which to spread reward (SD of Gaussian) (0..1)

%define radial basis
[~, ~, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(ntimesteps, nbasis, prop_spread);

%states for random generators are shared across functions to allow for repeatable draws
global rew_rng_state;

%initialize states for two repeatable random number generators using different seeds
rew_rng_seed=rngseeds(1); %used by RewFunction to draw probabilitic outcomes from contingency (not relevant for struct-based input)
explore_rng_seed=rngseeds(2); %Used by agents to choose between exploratory and exploitative choices. Only pertains to non-softmax agents.
grw_step_rng_seed=rngseeds(3); %used by GRW agent to sample random normal numbers for GRW step size
randrt_seed=rngseeds(4); %used by epsilon greedy agent with random uniform exploration
softmax_seed=rngseeds(5); %used by softmax agent for weighted sampling of softmax function using randsample

%populate free parameters for the requested agent/model
%note: Kalman filter does not have a free learning rate parameter.
if strcmpi(agent, 'fixedLR_softmax')
    %Learning: Bush-Mosteller fixed learning rate applied to PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: softmax function of value to select among alternatives
    beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)
    alpha = params(3);          %learning rate (0..1)
elseif strcmpi(agent, 'fixedLR_egreedy')
    %Learning: Bush-Mosteller fixed learning rate applied to PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: Epsilon-greedy with hard max of value for exploit and random uniform explore
    epsilon = params(2);        %proportion of exploratory choices (0..1)
    alpha = params(3);          %learning rate (0..1)
    
    rng(randrt_seed);
    randuniform_rts_explore = randi([min(tvec), max(tvec)], ntrials, 1); %vector of random uniform explore RTs
elseif strcmpi(agent, 'fixedLR_egreedy_grw')
    %Learning: Bush-Mosteller fixed learning rate applied to PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: Epsilon-greedy with hard max of value for exploit and GRW exploration
    epsilon = params(2);                %proportion of exploratory choices (0..1)
    alpha = params(3);                  %learning rate (0..1)
    sig_grw = params(4)*range(tvec);    %SD of GRW (0..1 interval proportion) rescaled wrt the observed interval
    
    %determine whether to strategic explore versus GRW (not currently used)
    %rng(exptype_rng_seed);
    %explore_type_rand=rand(ntrials,1);
    %exptype_rng_state=rng; %no need to save if this is not re-used.
    
    rng(grw_step_rng_seed); %seed for GRW step sizes
    grw_step=round(sig_grw*randn(ntrials,1)); %compute vector of step sizes for GRW
elseif strcmpi(agent, 'asymfixedLR_softmax')
    %Learning: Bush-Mosteller separate fixed learning rates for PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: softmax function of value to select among alternatives
    beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)
    alpha = params(3);          %learning rate for PE+ (0..1)
    rho = params(4);            %learning rate for PE- (0..1)
elseif strcmpi(agent, 'kalman_softmax')
    %Learning: kalman filter, gain initialized to 0.5 with Bayesian update (roughly exponential decay); uncertainty is tracked, but not used
    %Choice: softmax function of value to select among alternatives (no uncertainty)
    beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)  
elseif strcmpi(agent, 'kalman_processnoise')
    %Learning: kalman filter, gain initialized to 0.5 with Bayesian update (roughly exponential decay); uncertainty is tracked, but not used
    %Choice: softmax function of value to select among alternatives
    %Additional: PEs boost process noise (Q), effectively enhancing learning rate when surprising events occur.
    beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)  
    omega = params(3);          %scaling of process noise Q by PE.
elseif strcmpi(agent, 'kalman_sigmavolatility')
    %Learning: kalman filter, gain initialized to 0.5 with Bayesian update (roughly exponential decay); no representation of uncertainty
    %Choice: softmax function of value to select among alternatives
    %Additional: Prediction errors increase uncertainty (sigma) by a smooth function of PEs according to a second learning rate, phi.
    beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)  
    phi = params(3);            %scales additional noise added to sigma update according to smooth function of PEs.
    gamma = params(4);          %decay factor for volatility (0..1)
elseif strcmpi(agent, 'kalman_uv_logistic')
    %old logistic explore/exploit choice rule
    tradeoff = params(2); %point at which agent is indifferent between explore and exploit
    discrim = params(3); %discrimination/slope of logistic function for explore/exploit sampling (0.01..100)
elseif strcmpi(agent, 'kalman_uv_sum')
    beta = params(2);
    tau = params(3); % tau from Greek ???? -- value, price, cf ???? as trophys in Homer
elseif strcmpi(agent, 'fixedLR_kl_softmax')
    beta = params(2);
    alpha = params(3);   %learning rate for value
    kappa = params(4);   %PE+ tilt scaling parameter
    lambda = params(5);  %PE- tilt scaling parameter 
elseif strcmpi(agent, 'kalman_kl_softmax')
    beta = params(2);
    kappa = params(3);
    lambda = params(4);
elseif strcmpi(agent, 'kalman_processnoise_kl')
    beta = params(2);
    omega = params(3);  %scaling of process noise by PE.
    kappa = params(4);  %PE+ tilt scaling parameter
    lambda = params(5); %PE- tilt scaling parameter
elseif strcmpi(agent, 'kalman_uv_sum_kl')
    beta = params(2);
    tau = params(3);
    kappa = params(4);
    lambda = params(5);
end

if prop_spread < 0 || prop_spread > 1, error('prop_spread outside of bounds'); end

%cond can be a character string (e.g., DEV) in which case the agent draws from the reward function on each trial.
%it can also be a struct, in which case this defines the lookup table for sequential samples from a contingency of interest.
usestruct=0;
if isstruct(cond)
    cstruct = cond;
    cond = cstruct.name;
    usestruct=1;
end

%add Gaussian noise with sigma = 1% of the range of the time interval to rt_explore
%prop_expnoise=.01;
%sig_expnoise=prop_expnoise*range(tvec);

%% initialize RTs chosen by agent as a nan vector;
rts = nan(1,ntrials);

%Have initial rts be constant for cost function
rts(1) = ceil(.5*ntimesteps);

%setup matrices for tracking learning
% i = trial
% j = basis function
% t = time step within trial, in centiseconds (1-500, representing 0-5 seconds)
delta_ij =      zeros(ntrials, nbasis);     % prediction error assigned to each microstimulus
e_ij =          zeros(ntrials, nbasis);     % eligibility traces for each microstimulus in relation to RT (US)
v_jt =          zeros(nbasis, ntimesteps);  % value by microstimulus (rows for every microstimulus and columns for time points within trial)
v_it =          zeros(ntrials, ntimesteps); % history of value function by trial
vfinal_it =     zeros(ntrials, ntimesteps); % history of value function by trial after manipulations (u+v etc.)
rew_i =         nan(1, ntrials);            % actual reward for each trial
ev_i =          nan(1, ntrials);            % expected value of choice for each trial (used for cost function)
u_jt =          zeros(nbasis, ntimesteps);  % uncertainty of each basis for each timestep
u_it =          zeros(ntrials,ntimesteps);  % history of uncertainty by trial at each timestep
uv_it =         zeros(ntrials,ntimesteps);  % history of UV (value incorporating exploration bonus) by trial at each timestep
z_i =           nan(1, ntrials);            % smooth volatility function in kalman_sigmavolatility model (not wrt basis function)
%delta_func =    zeros(1,ntrials);
%v_max_to_plot=  zeros(1,ntrials);
p_choice    =   nan(ntrials,ntimesteps);    % model-predicted probabilities of choosing each possible rt according to softmax

%kalman setup
mu_ij =         nan(ntrials, nbasis);       %means of Gaussians for Kalman
k_ij =          zeros(ntrials, nbasis);     %Kalman gain (learning rate)
sigma_ij =      zeros(ntrials, nbasis);     %Standard deviations of Gaussians (uncertainty)
Q_ij =          zeros(ntrials, nbasis);     %Process noise (i.e., noise/change in the system dynamics with time)

mu_ij(1,:) =    0; %expected reward on first trial is initialized to 0 for all Gaussians.
z_i(1) = 0;        %no initial expected volatility. 

%despair = zeros(1,ntrials);                 %despair/volatility to detect reversals

%noise in the reward signal: sigma_rew. In the Frank model, the squared SD of the Reward vector represents the noise in
%the reward signal, which is part of the Kalman gain. This provides a non-arbitrary initialization for sigma_rew, such
%that the variability in returns is known up front... This is implausible from an R-L perspective, but not a bad idea
%for getting a reasonable estimate of variability in returns. It would be interesting (plausible) to give a unique
%estimate to each basis based on its temporal receptive field, but hard to decide this up front. But it does seem like
%this could give local sensitivity to volatility (e.g., low- versus high-probability windows in time). For now, I'll
%just use the variance of the returns (ala Frank). But for the generative agent, this is not known up front -- so sample
%from the chosen contingency for each timestep as a guess.

if usestruct
    %taking std over all timesteps and possible draws here. This is in contrast to approach below where you get one
    %sample from the contingency as drawn by RewFunction. This is much faster than arrayfun below.
    sigma_noise = repmat(std(cstruct.lookup(:))^2, 1, nbasis); 
else
    rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
    rew_rng_state=rng;
    sigma_noise = repmat(std(arrayfun(@(x) RewFunction(x*10, cond, 0), tvec))^2, 1, nbasis);
end

%Define indifference point between explore and exploit (p = 0.5) as proportion reduction in variance from initial value
u_threshold = (1-tradeoff) * sigma_noise(1);

%As in Frank, initialize estimate of std of each Gaussian to the noise of returns on a sample of the whole contingency.
%This leads to an effective learning rate of 0.5 since k = sigma_ij / sigma_ij + sigma_noise
sigma_ij(1,:) = sigma_noise;

%fprintf('running agent with sigs: %.3f, tradeoff: %.3f and rngseeds: %s \n', sig_spread, tradeoff, num2str(rngseeds));

%rng calls are slow according to profiler
%pull a vector of random numbers in advance of the trial loop and pull ith element

%generate vector of random uniform numbers to define random/exploratory choice probabilities for each trial
%this is only used by agents that do not use the softmax choice function (epsilon greedy and uv logistic)
if ismember(agent, {'kalman_uv_logistic', 'fixedLR_egreedy', 'fixedLR_egreedy_grw'})
    rng(explore_rng_seed); %draw from explore/exploit rng
    choice_rand=rand(ntrials,1);
end

%setup random number generator for softmax function if relevant
if ismember(agent, {'fixedLR_softmax', 'asymfixedLR_softmax', 'kalman_softmax', 'kalman_processnoise', ...
        'kalman_sigmavolatility', 'kalman_uv_sum', 'fixedLR_kl_softmax', ...
        'kalman_kl_softmax', 'kalman_processnoise_kl', 'kalman_uv_sum_kl'})
    
    %setup a random number stream to be used by the softmax choice rule to make choices consistent during optimization
    softmax_stream = RandStream('mt19937ar','Seed',softmax_seed);
end 

%Set up to run multiple runs for multiple ntrials
for i = 1:ntrials    
    % get symmetric eligibility traces for each basis function (temporal generalization)
    % generate a truncated Gaussian basis function centered at the RT and with sigma equal to the free parameter.
    
    %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
    elig = gaussmf(tvec, [sig_spread, rts(i)]);
    
    %compute sum of area under the curve of the gaussian function
    auc=sum(elig);
    
    %divide gaussian update function by its sum so that AUC=1.0, then rescale to have AUC of a non-truncated basis
    %this ensures that eligibility is 0-1.0 for non-truncated update function, and can exceed 1.0 at the edge.
    %note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
    %will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
    elig=elig/auc*refspread;
    
    %truncated gaussian eligibility
    %figure(7); plot(tvec, elig);
    
    %compute the intersection of the Gaussian spread function with the truncated Gaussian basis.
    %this is essentially summing the area under the curve of each truncated RBF weighted by the truncated
    %Gaussian spread function.
    e_ij(i,:) = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 2);
    
    %%NB: broken for now
%     if reversal==1 && i==ntrials/2+1
%         %Optimal param reversal hack
%         vperm_run = m.vperm_run; %Currently this only exsists for the reversal runs
%         if strcmp(m.name, 'IEV')
%             m=mDEV; %if it is IEV after x trials switch
%             m.lookup = m.lookup(:,vperm_run);
%             cond = m.name;
%         else
%             m=mIEV; %else it is DEV after x traisl switch to IEV
%             m.lookup = m.lookup(:,vperm_run);
%             cond = m.name;
%         end        
%     end

    if usestruct
        [rew_i(i), ev_i(i), cstruct] = getNextRew(rts(i), cstruct);
        
        %The code below is slow due to rand/rng usage and computing functions, rather than struct pre-compute approach.
        %[~, ev_i(i)] = RewFunction(rts(i).*10, cond, 0); %multiply by 10 because underlying functions range 0-5000ms
        %don't use the reward function rng seed if already passed in master structure (since we're just looking up EV).
    else
        [rew_i(i), ev_i(i)] = RewFunction(rts(i).*10, cond); %multiply by 10 because underlying functions range 0-5000ms
    end
    
    %1) compute prediction error, scaled by eligibility trace
    delta_ij(i,:) = e_ij(i,:).*(rew_i(i) - mu_ij(i,:));

    %Variants of learning rule
    if ismember(agent, {'fixedLR_softmax', 'fixedLR_egreedy', 'fixedLR_egreedy_grw', 'fixedLR_kl_softmax'})
        mu_ij(i+1,:) = mu_ij(i,:) + alpha.*delta_ij(i,:);
    elseif strcmpi(agent, 'asymfixedLR_softmax')
        %need to avoid use of mean function for speed in optimization... would max work?
        if (max(delta_ij(i,:))) > 0
            mu_ij(i+1,:) = mu_ij(i,:) + alpha.*delta_ij(i,:); %learn from PE+
        else
            mu_ij(i+1,:) = mu_ij(i,:) + rho.*delta_ij(i,:); %learn from PE-
        end
    else
        %Kalman variants of learning and uncertainty
        
        %Changing Kalman variance a posteriori uses the elig*gain approach: [1 - k(ij)*elig(ij)]*sigma(ij)
        %this would only allow a 1.0 update*kalman gain for basis functions solidly in the window and a decay in diminishing
        %variance as the basis deviates from the timing of the obtained reward.
        
        %MH 8Sep2015: At the moment, we assume zero process noise in the estimated posterior error covariances, sigma_ij.
        %To model dynamic/changing systems, try dynamically enhance learning rates by scaling process noise by PE.
        Q_ij(i,:) = omega.*abs(delta_ij(i,:)); %use abs of PE so that any large surprise enhances effective gain.
        
        %Compute the Kalman gains for the current trial (potentially adding process noise)
        k_ij(i,:) = (sigma_ij(i,:) + Q_ij(i,:))./(sigma_ij(i,:) + Q_ij(i,:) + sigma_noise);       
        
        %Update posterior variances on the basis of Kalman gains
        sigma_ij(i+1,:) = (1 - e_ij(i,:).*k_ij(i,:)).*(sigma_ij(i,:) + z_i(i));
        
        %Update reward expectation. AD: Would it be better for the delta to be the difference between the reward
        %and the point value estimate at the RT(i)?
        mu_ij(i+1,:) = mu_ij(i,:) + k_ij(i,:).*delta_ij(i,:);
        
        %Track smooth estimate of volatility according to unsigned PE history
        z_i(i+1) = gamma.*z_i(i) + phi.*abs(sum(delta_ij(i,:)));
        
        %Uncertainty is a function of Kalman uncertainties.
        u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat;        
        u_func = sum(u_jt); %vector of uncertainties by timestep
        u_it(i+1,:) = u_func;
    end
        
    %compute summed/evaluated value function across all timesteps
    v_jt=mu_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
    v_func = sum(v_jt); %subjective value by timestep as a sum of all basis functions
    
    %% How does the agent detect a reversal?  E.g. by seeing a lot of PE-.  These, in turn, inflate a
    %  'despair' parameter (analogous to Behrens's volatility [2005]).
    % first get the PE-(i)
    % removing despair computation for now since it is not used and mean call is relatively slow
    %mean_delta(i) = mean(delta_ij(i,:));
    %mean_delta_minus(i) = mean_delta(i)*(mean_delta(i)<0);
    %despair(i+1) = .95*despair(i) - mean_delta_minus(i);
    
    %Prediction errors
    %delta_func(i)=sum(delta_ij(i,:));
    %v_max_to_plot(i) = max(v_func)*50;
    %k_top_plot(i) = sum(k_ij(i,:));
    
    if i == ntrials, break; end %do not compute i+1 choice on the final trial (invalid indexing problem)
    
    %compute hard max of value function alone (used by kalman_uv_logistic, fixedLR_egreedy, and fixedLR_egreedy_grw)
    if sum(v_func) == 0
        rt_exploit = ceil(.5*ntimesteps); %default to mid-point of time domain
    else
        rt_exploit = find(v_func==max(v_func), 1); %only take the first max if there is are two identical peaks.
        if rt_exploit > max(tvec), rt_exploit = max(tvec); end
    end
    
    %% Choice rule
    if strcmpi(agent, 'kalman_uv_logistic')
        %compared to other models that use a curve over which to choose (either by softmax or egreedy selection),
        %kalman_uv_logistic computes explore and exploit choices and chooses according to a logistic.
        u = sum(u_func)/length(u_func);

        if u == 0
            rt_explore = ceil(.5*ntimesteps);
        else
            rt_explore = find(u_func==max(u_func), 1); %return position of first max (and add gaussian noise?)
        end

        sigmoid = 1/(1+exp(-discrim.*(u - u_threshold))); %Rasch model with tradeoff as difficulty (location) parameter
        
        if choice_rand(i) < sigmoid
            %explore according to hardmax u
            rts(i+1) = rt_explore;
        else
            rts(i+1) = rt_exploit;
        end
        
        v_final = v_func; %no alterations of value function for logistic
    else
        %compute kappa/lambda discount vector to be added to decision function below
        if ismember(agent, {'fixedLR_kl_softmax', 'kalman_kl_softmax', 'kalman_processnoise_kl', 'kalman_uv_sum_kl'})
            %compute a "tilt" vector that linearly scales with timepoint according to PE size * parameter
            %use different parameters, kappa and lambda, to handle PE+ and PE-, respectively
            if max(delta_ij(i,:)) > 0
                discount = kappa*(max(delta_ij(i,:)))*tvec; %% also add valence-dependent parameters: kappa for PE+, lambda for PE-
            else
                discount = lambda*(min(delta_ij(i,:)))*tvec;
            end
        end
        
        %compute final value function to use for choice
        if ismember(agent, {'fixedLR_softmax', 'fixedLR_egreedy', 'fixedLR_egreedy_grw', ...
                'asymfixedLR_softmax', 'kalman_softmax', 'kalman_processnoise', 'kalman_sigmavolatility'})
            v_final = v_func; % just use value curve for choice
        elseif strcmpi(agent, 'kalman_uv_sum')
            uv_func=tau*v_func + (1-tau)*u_func; %mix together value and uncertainty according to tau
            v_final = uv_func;
            uv_it(i+1,:) = uv_func;
        elseif ismember(agent, {'kalman_kl_softmax', 'fixedLR_kl_softmax', 'kalman_processnoise_kl'})
            v_final= v_func + discount;
            %v_it_undisc(i+1,:) = v_func; %not implemented here yet
            
            %not storing rt_exploit_disc yet (see line 372 of skeptic_fitsubject_all_models.m)
        elseif strcmpi(agent, 'kalman_uv_sum_kl')
            v_disc = v_func + discount;
            %v_it_undisc(i+1,:) = v_func;
            uv_func=tau*v_disc + (1-tau)*u_func;
            v_final = uv_func;
        end
        
        %compute choice rule according to agent
        if ismember(agent, {'fixedLR_egreedy', 'fixedLR_egreedy_grw'})
            if choice_rand(i) < epsilon %explore
                if strcmpi(agent, 'fixedLR_egreedy')
                    rts(i+1) = randuniform_rts_explore(i);
                elseif strcmpi(agent, 'fixedLR_egreedy_grw')
                    rt_grw = rts(i) + grw_step(i);
                
                    %N.B.: Need to have more reasonable GRW near the edge such that it doesn't just oversample min/max
                    %e.g., perhaps reflect the GRW if rt(t-1) was already very close to edge and GRW samples in that direction again.
                    if rt_grw > max(tvec), rt_grw = max(tvec);
                    elseif rt_grw < min(tvec), rt_grw = min(tvec); end
                    rts(i+1) = rt_grw;
                end
            else
                rts(i+1) = rt_exploit; %need to unify this with kalman_uv_logistic above...
            end
        else
            %ismember(agent, {'fixedLR_softmax', 'asymfixedLR_softmax', 'kalman_softmax'})
            %NB: all other models use a softmax choice rule over the v_final curve.
            p_choice(i,:) = (exp((v_final-max(v_final))/beta))/(sum(exp((v_final-max(v_final))/beta))); %Divide by temperature
            %if (all(v_final==0)), v_final=rand(1, length(v_final)).*1e-6; end; %need small non-zero values to unstick softmax on first trial
            rts(i+1) = randsample(softmax_stream, tvec, 1, true, p_choice(i,:));
        end
        
    end
    
    %populate v_it for tracking final value function
    vfinal_it(i+1,:) = v_final; %store choice function for return according to model

    %Output the rt explore and exploit from choice rule per trial
    %fprintf('trial: %d rt_exploit: %.2f rt_explore: %.2f\n', i, rt_exploit, rt_explore);
        
    verbose=0;
    if verbose == 1
        fprintf('Trial: %d, Rew(i): %.2f, Rt(i): %.2f\n', i, rew_i(i), rts(i));
        %fprintf('w_i,k:    '); fprintf('%.2f ', mu_ij(i,:)); fprintf('\n');
        %fprintf('delta_ij:   '); fprintf('%.2f ', delta_ij(i,:)); fprintf('\n');
        %fprintf('w_i+1,k:  '); fprintf('%.2f ', mu_ij(i+1,:)); fprintf('\n');
        fprintf('\n');        
    end
    
end

%Cost function
cost = -sum(ev_i);

ret.agent=agent;
ret.nbasis = nbasis;
ret.ntimesteps = ntimesteps;
ret.sigma_noise=sigma_noise;
ret.rts=rts;
ret.mu_ij = mu_ij;
ret.sigma_ij = sigma_ij;
ret.k_ij = k_ij;
ret.delta_ij = delta_ij;
ret.Q_ij = Q_ij;
ret.e_ij = e_ij;
ret.v_it = v_it;
ret.vfinal_it = vfinal_it;
ret.u_it = u_it;
ret.uv_it = uv_it;
ret.rew_i = rew_i;
ret.ev_i = ev_i;
ret.z_i = z_i;
ret.p_choice=p_choice;