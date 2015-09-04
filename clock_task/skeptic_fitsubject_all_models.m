%Kalman filter extension of temporal basis operator: each basis function tracks the mean value and uncertainty (sd) as a
%Gaussian.
%Authors: Michael Hallquist, Alex Dombrovski, & Jonathan Wilson
%Date last modified: 4/21/2015
%Matlab Version: R2012a
%fit subject behavior using skeptic (Strategic Kalman filter Exploration/Exploitation of Temporal Instrumental
%Contingencies)

%% TODO
% 1) downweight costs at beginning of the trial prior to when people are likely able to respond (< 250ms)
%       - This may just mean setting a minimum acceptable model-predicted RT (censor)


%%
function [cost,ret,mov] = skeptic_fitsubject_all_models(params, rt_obs, rew_obs, rngseeds, nbasis, ntimesteps, trial_plots, minrt, maxrt, modelname)
%params is the vector of parameters used to fit
%   params(2): prop_spread -- width of temporal generalization Gaussian as a proportion of the time interval (0..1)
%   params(3): k -- proportion of strategic versus stochastic exploration
%   params(4): s_grw -- width of Gaussian random walk (GRW) for stochastic exploration
%cond is the character string of the reward contingency
%ntrials is the number of trials to run
%nbasis is the number of radial basis functions used to estimate value and uncertainty
%ntimesteps is the number of time bins used for obtaining estimates of time functions for plotting etc.

ret=[];
ntrials=length(rew_obs);

if nargin < 4
    rngseeds=[98 83 66 10];
end

if nargin < 5, nbasis = 24; end
if nargin < 6, ntimesteps=400; end
if nargin < 7, trial_plots = 1; end
if nargin < 8, minrt=25; end %minimum 250ms RT
if nargin < 9, maxrt=400; end
if nargin < 10, modelname='value_softmax'; end

%note: Kalman filter does not have a free learning rate parameter.
if strcmpi(modelname, 'value_softmax')
    %% free parameters:
    %    1) prop_spread: width of temporal generalization function as a proportion of interval. 0--1
    %    2) beta: inverse temperature of softmax. 0--infinity
    prop_spread = params(1);
    beta = params(2);
    
    
elseif strcmpi(modelname, 'uv')
    %% free params
    prop_spread = params(1);
    tau = params(2); % tau from Greek ???? -- value, price, cf ???? as trophys in Homer
    beta = params(3);
    
elseif strcmpi(modelname, 'v_discounted')
    prop_spread = params(1);
    beta = params(2);
    kappa = params(3);
elseif strcmpi(modelname, 'uv_discounted')
    prop_spread = params(1);
    beta = params(2);
    kappa = params(3);
    tau = params(4);
end

if prop_spread < 0 || prop_spread > 1, error('prop_spread outside of bounds'); end


%treats time steps and max time as synonymous
fprintf('Downsampling rt_obs by a factor of 10 to to 0..%d\n', ntimesteps);
rt_obs = round(rt_obs/10);

%subtract off minimum rt so that subject's min RT or 25 (250ms) is now treated as 0
%also need to shift the basis by the corresponding amount so that
%ntimesteps = 375 if min RT = 250ms; also truncate basis on the right (late)
%end to reflect subject's max RT
rt_obs = rt_obs - minrt+1;
%occasional subjects are very fast, so
rt_obs(rt_obs<1) = 1;
ntimesteps = maxrt - minrt + 1;

%states for random generators are shared across functions to allow for repeatable draws
global rew_rng_state explore_rng_state;

%initialize states for two repeatable random number generators using different seeds
rew_rng_seed=rngseeds(1);
explore_rng_seed=rngseeds(2);
grw_step_rng_seed=rngseeds(3);
exptype_rng_seed=rngseeds(4);

rng(rew_rng_seed);
rew_rng_state=rng;
rng(explore_rng_seed);
explore_rng_state=rng;
rng(grw_step_rng_seed);
grw_step_rng_state=rng;
rng(exptype_rng_seed);
exptype_rng_seed=rng;

%initialize movie storage
mov=repmat(struct('cdata', [], 'colormap', []), ntrials,1);

%Initialize time step vector and allocate for memory
tvec=1:ntimesteps;
sig_spread=prop_spread*range(tvec); %determine SD of spread function

%% GRW not utilized in subject fitting
% %rescale s_grw wrt the interval (not as a proportion)
% s_grw=s_grw*range(tvec); %determine SD of spread function

% %add Gaussian noise with sigma = 1% of the range of the time interval to rt_explore
% prop_expnoise=.01;
% sig_expnoise=prop_expnoise*range(tvec);

%setup centers (means) and sds of basis functions
%based on testing in fix_rbf_basis.m, place the lowest center 12.5% below the first timestep and the last
%center 12.5% above last timestep. SD should be calculated to give a Cohen's d of 1.52 between
%basis functions (~45% distribution overlap).

%margin_offset=0;
margin_offset = (max(tvec) - min(tvec))*.125; % 12.5% offset

%define lowest and highest centers
tmin = min(tvec) - margin_offset; tmax=max(tvec) + margin_offset;
c=tmin:(tmax-tmin)/(nbasis-1):tmax;

sig = (c(2) - c(1))/1.52; %cohen's d of 1.52 between basis functions

%setup matrices for tracking learning
% i = trial
% j = basis function
% t = time step within trial, in centiseconds (1-500, representing 0-5 seconds)
delta_ij =      zeros(ntrials, nbasis);     % prediction error assigned to each microstimulus
e_ij =          zeros(ntrials, nbasis);     % eligibility traces for each microstimulus in relation to RT (US)
v_jt =          zeros(nbasis, ntimesteps);  % value by microstimulus (rows for every microstimulus and columns for time points within trial)
v_it =          zeros(ntrials, ntimesteps); % history of value function by trial
ev_obs_i =      nan(1, ntrials);            % expected value of observed choices for each trial (used for comparing subject and model)
ev_pred_i =     nan(1, ntrials);            % expected value of predicted choices
u_jt =          zeros(nbasis, ntimesteps);  % uncertainty of each basis for each timestep
u_it =          zeros(ntrials,ntimesteps);  % history of uncertainty by trial at each timestep
d_i =           nan(1, ntrials);            % euclidean distance (in U and Q space) between obs and pred RTs
precision_i =   nan(1, ntrials);            % model-predicted certainty about choice type (explore/exploit)... Used for downweighting costs for uncertain choices
p_choice    =   nan(ntrials,ntimesteps);    % model-predicted probabilities of choosing each possible rt according to softmax
p_chosen    =   nan(1, ntrials);            % model-predicted probabilities of the chosen option
uv_it =         zeros(ntrials,ntimesteps);  % history of UV (value incorporating exploration bonus) by trial at each timestep
v_it_undisc =   zeros(ntrials, ntimesteps); % history of value function by trial, undiscounted for opportunity cost (for v_discounted model)


%kalman setup
mu_ij =         nan(ntrials, nbasis);       %means of Gaussians for Kalman
k_ij =          zeros(ntrials, nbasis);     %Kalman gain (learning rate)
sigma_ij =      zeros(ntrials, nbasis);     %Standard deviations of Gaussians (uncertainty)

mu_ij(1,:) =    0; %expected reward on first trial is initialized to 0 for all Gaussians.

rt_pred_i =          nan(1, ntrials);            % predicted reaction times
rts_pred_explore =   nan(1, ntrials);            % predicted exploratory reaction times
rts_pred_exploit =   nan(1, ntrials);            % predicted exploitative reaction times
exploit_trials =     [];                         %vector of trials where model predicts exploitation
explore_trials =     [];                         %vector of trials where model predicts exploration
p_explore_i =        nan(1, ntrials);            %model-predicted probability of exploration
rts_uv_pred =         nan(1, ntrials);            % predicted UV reaction times
rts_pred_exploit_disc =   nan(1, ntrials);            % predicted discounted exploitative reaction times
rts_pred_uv_disc =    nan(1, ntrials);            % predicted discounted exploitative reaction times

rt_pred_i(1) = rt_obs(1); %first choice is exogenous to model
rts_pred_explore(1) = rt_obs(1); %first choice is exogenous to model
rts_pred_exploit(1) = rt_obs(1); %first choice is exogenous to model
rts_uv_pred(1) = rt_obs(1); %first choice is exogenous to model

%noise in the reward signal: sigma_rew. In the Frank model, the squared SD of the Reward vector represents the noise in
%the reward signal, which is part of the Kalman gain. This provides a non-arbitrary initialization for sigma_rew, such
%that the variability in returns is known up front... This is implausible from an R-L perspective, but not a bad idea
%for getting a reasonable estimate of variability in returns. It would be interesting (plausible) to give a unique
%estimate to each basis based on its temporal receptive field, but hard to decide this up front. But it does seem like
%this could give local sensitivity to volatility (e.g., low- versus high-probability windows in time). For now, I'll
%just use the variance of the returns (ala Frank). But for the generative agent, this is not known up front -- so sample
%from the chosen contingency for each timestep as a guess.

%sigma_noise = repmat(std(arrayfun(@(x) RewFunction(x*10, cond, 0), tvec))^2, 1, nbasis);

%use observed volatility of choices (gives model some unfair insight into future... but need a basis for
%expected volatility/process noise)
sigma_noise = repmat(var(rew_obs), 1, nbasis);

%As in Frank, initialize estimate of std of each Gaussian to the noise of returns on a sample of the whole contingency.
%This leads to an effective learning rate of 0.5 since k = sigma_ij / sigma_ij + sigma_noise
sigma_ij(1,:) = sigma_noise;

%so, Kalman Gaussians should be updated by the obtained reward (alternative: PE), spread in time as usual, which forms a
%temporal eligibility trace... But how do we make sure that we don't get big PEs by gain*(Rew - mean)? Should be a
%similar problem to the existing guy. In the current implementation we use the alpha*elig*Reward - Expected), such that
%we only effectively update value representation for large values of the eligibility trace. Should be parallel here:
%elig*gain*(Rew - mean).

%construct radial basis matrix using Gaussians
gaussmat = zeros(nbasis,ntimesteps);

for j = 1:nbasis
    gaussmat(j,:) = gaussmf(tvec,[sig c(j)]);
end

%version of gaussian where each function has AUC = 1.0 (PDF representation)
maxauc_all=max(sum(gaussmat, 2));
gaussmat_pdf=gaussmat./maxauc_all;

%normalize gauss functions to each have AUC = 1.0 within observed time interval
%this is essentially a truncated Gaussian basis such that AUC = 1.0 for basis functions within interval
maxauc_each=sum(gaussmat,2)*ones(1,length(tvec)); %outer product of vectors to allow for col-wise division below
gaussmat_trunc=gaussmat./maxauc_each;

rbf_plots = 0;
if trial_plots == 1 && rbf_plots == 1
    figure(20); plot(tvec,gaussmat); title('Regular RBF');
    figure(21); plot(tvec,gaussmat_trunc); title('Truncated RBF');
end

%fprintf('updating value by alpha: %.4f\n', alpha);
%fprintf('updating value by epsilon: %.4f with rngseeds: %s \n', epsilon, num2str(rngseeds));
%fprintf('running agent with sigs: %.3f, epsilon: %.3f and rngseeds: %s \n', sig_spread, epsilon, num2str(rngseeds));

%determine the AUC of a non-truncated eligilibity function
refspread = sum(gaussmf(min(tvec)-range(tvec):max(tvec)+range(tvec), [sig_spread, median(tvec)]));

%objective expected value for this function
%ev=[];
%for val = 1:length(tvec)
%    [~,ev(val)] = RewFunction(tvec(val).*10, cond);
%end

%figure(6); plot(tvec, ev);
%title('Expected value of contingency');

%Set up to run multiple runs for multiple ntrials
for i = 1:ntrials
    
    % get symmetric eligibility traces for each basis function (temporal generalization)
    % generate a truncated Gaussian basis function centered at the RT and with sigma equal to the free parameter.
    
    %compute gaussian spread function with mu = rt_obs(i) and sigma based on free param prop_spread
    elig = gaussmf(tvec, [sig_spread, rt_obs(i)]);
    
    %compute sum of area under the curve of the gaussian function
    auc=sum(elig);
    
    %divide gaussian update function by its sum so that AUC=1.0, then rescale to have AUC of a non-truncated basis
    %this ensures that eligibility is 0-1.0 for non-truncated update function, and can exceed 1.0 at the edge.
    %note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
    %will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
    elig=elig/auc*refspread;
    
    %truncated gaussian eligibility
    %figure(7); plot(tvec, elig);
    
    %compute the product of the Gaussian spread function with the truncated Gaussian basis.
    e_ij(i,:) = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 2);
    
    %Changing Kalman variance a posteriori should also use the elig*gain approach: [1 - k(ij)*elig(ij)]*sigma(ij)
    %this would only allow a 1.0 update*kalman gain for basis functions solidly in the window and a decay in diminishing
    %variance as the basis deviates from the timing of the obtained reward.
    
    %1) compute the Kalman gains for the current trial
    k_ij(i,:) = sigma_ij(i,:)./(sigma_ij(i,:) + sigma_noise);
    
    %2) update posterior variances on the basis of Kalman gains
    sigma_ij(i+1,:) = (1 - e_ij(i,:).*k_ij(i,:)).*sigma_ij(i,:);
    
    %3) update reward expectation
    delta_ij(i,:) = e_ij(i,:).*(rew_obs(i) - mu_ij(i,:));
    mu_ij(i+1,:) = mu_ij(i,:) + k_ij(i,:).*delta_ij(i,:);
    
    v_jt=mu_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
    
    %subjective value by timestep as a sum of all basis functions
    v_func = sum(v_jt);
    
    %uncertainty is now a function of Kalman uncertainties.
    u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat;
    u_func = sum(u_jt);
    
    v_it(i+1,:) = v_func;
    u_it(i+1,:) = u_func;
    
    %compute rt_exploit
    if sum(v_func) == 0
        rt_exploit = rt_obs(1); %feed RT exploit the first observed RT
    else
        rt_exploit = find(v_func==max(v_func));
        
        if rt_exploit > max(tvec)
            rt_exploit = max(tvec);
        elseif rt_exploit < minrt
            rt_exploit = minrt; %do not allow choice below earliest possible response
        end
    end
    %% just for plotting
    rts_pred_exploit(i+1) = rt_exploit;
    
    %% compute final value function used for decision
    if strcmpi(modelname, 'value_softmax')
        %rt_explore = rt_exploit + grw_step;
        v_final = v_func;
        
    elseif strcmpi(modelname, 'uv')
        uv_func=tau*v_func + (1-tau)*u_func;
        v_final = uv_func;
        uv_it(i+1,:) = uv_func;
    elseif strcmpi(modelname, 'v_discounted')
        vec = 1:ntimesteps;
        discount = kappa*(mean(delta_ij(i,:)))*vec;
        v_final = v_func + discount;
        v_it_undisc(i+1,:) = v_func;
        v_it(i+1,:) = v_final;
        %% RTexploit for plotting:
        if sum(v_func) == 0
            rt_exploit_disc = rt_obs(1); %feed RT exploit the first observed RT
        else
            rt_exploit_disc = find(v_final==max(v_final));
            
            if rt_exploit_disc > max(tvec)
                rt_exploit_disc = max(tvec);
            elseif rt_exploit_disc < minrt
                rt_exploit_disc = minrt; %do not allow choice below earliest possible response
            end
        end
        rts_pred_exploit_disc(i+1) = rt_exploit_disc;
    elseif strcmpi(modelname, 'uv_discounted')
        vec = 1:ntimesteps;
        discount = kappa*(mean(delta_ij(i,:)))*vec;
        v_final = v_func + discount;
        v_it_undisc(i+1,:) = v_func;
        v_it(i+1,:) = v_final;
        uv_func=tau*v_final + (1-tau)*u_func;
        v_final = uv_func;
        uv_it(i+1,:) = uv_func;
        %% RTexploit for plotting:
        if sum(v_func) == 0
            rt_exploit_disc = rt_obs(1); %feed RT exploit the first observed RT
        else
            rt_exploit_disc = find(v_it(i+1,:)==max(v_it(i+1,:)));
            
            if rt_exploit_disc > max(tvec)
                rt_exploit_disc = max(tvec);
            elseif rt_exploit_disc < minrt
                rt_exploit_disc = minrt; %do not allow choice below earliest possible response
            end
        end
        rts_pred_exploit_disc(i+1) = rt_exploit_disc;
        %% RT_UV for plotting
        if sum(v_func) == 0
            rt_uv_disc = rt_obs(1); %feed RT exploit the first observed RT
        else
            rt_uv_disc = find(uv_func == max(uv_func));
            if rt_uv_disc > max(tvec)
                rt_uv_disc = max(tvec);
            elseif rt_uv_disc < minrt
                rt_uv_disc = minrt; %do not allow choice below earliest possible response
            end
        end
    end
        rts_pred_uv_disc(i+1) = rt_uv_disc;
        p_choice(i,:) = (exp((v_final-max(v_final)))/beta)/(sum(exp((v_final-max(v_final)))/beta)); %Divide by temperature
        
        p_chosen(i) = p_choice(i,rt_obs(i));
        
        
    end
    
    cost = -log(prod(p_chosen)); %maximize choice probabilities
    ret.mu_ij = mu_ij;
    ret.sigma_ij = sigma_ij;
    ret.k_ij = k_ij;
    ret.delta_ij = delta_ij;
    ret.e_ij = e_ij;
    ret.v_it = v_it;
    ret.u_it = u_it;
    ret.rew_obs = rew_obs;
    ret.ev_obs_i = ev_obs_i;
    ret.ev_pred_i = ev_pred_i;
    ret.rt_obs = rt_obs;
    ret.rt_pred_i = rt_pred_i;
    ret.rts_pred_explore = rts_pred_explore;
    ret.rts_pred_exploit = rts_pred_exploit;
    ret.explore_trials = explore_trials;
    ret.exploit_trials = exploit_trials;
    ret.d_i = d_i;
    ret.p_explore_i = p_explore_i;
    ret.p_choice = p_choice;
    ret.rt_uv_pred = rts_uv_pred;
    ret.p_chosen = p_chosen;
    ret.uv = uv_it;
    ret.v_it_undisc = v_it_undisc;
    ret.rts_pred_exploit_disc = rts_pred_exploit_disc;
    ret.rts_pred_uv_disc = rts_pred_uv_disc;
    
    if trial_plots
        %         figure(1); clf;
        %         subplot(5,2,1)
        %         %plot(tvec,v_func);
        %         scatter(rt_pred_i(1:i),rew_obs(1:i)); axis([1 500 0 350]);
        %         hold on;
        %         plot(rt_pred_i(i),rew_obs(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
        %         hold off;
        %         subplot(5,2,2)
        %         plot(tvec,v_func); xlim([-1 ntimesteps+1]);
        %         ylabel('value')
        %         subplot(5,2,3)
        %
        % %         bar(c, mu_ij(i,:));
        % %         ylabel('basis function heights');
        %         plot(tvec,v_jt);
        %         ylabel('temporal basis function')
        % %         title(sprintf('trial # = %i', h)); %
        %                 xlabel('time(ms)')
        %                 ylabel('reward value')
        %
        %         subplot(5,2,4)
        %         plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
        %         xlabel('time (centiseconds)')
        %         ylabel('uncertainty')
        %
        %         subplot(5,2,5)
        %         barh(sigmoid);
        %         xlabel('p(explore)'); axis([-.1 1.1 0 2]);
        %         subplot(5,2,6)
        %         %barh(alpha), axis([0 .2 0 2]);
        %         %xlabel('learning rate')
        %         subplot(5,2,7)
        %         barh(sig_spread), axis([0.0 0.5 0 2]);
        %         xlabel('decay')
        %         subplot(5,2,8)
        %         barh(epsilon), axis([-0.5 0 0 2]);
        %         xlabel('strategic exploration')
        %
        %         subplot(5,2,9)
        %         barh(u) %, axis([0 1000 0 2]);
        %         xlabel('mean uncertainty')
        %         %         pause(0.1);
        %         subplot(5,2,10)
        %         plot(1:ntrials,rt_pred_i, 'k');
        %         ylabel('rt by trial'); axis([1 ntrials -5 505]);
        
        figure(1); clf;
        set(gca,'FontSize',24);
        %         subplot(3,2,1);
        %         title('Choice history');
        %         %plot(tvec,v_func);
        %         scatter(rt_obs(1:i),rew_obs(1:i)); axis([1 ntimesteps 0 350]);
        %         hold on;
        %         plot(rt_obs(i),rew_obs(i),'r*','MarkerSize',20);  axis([1 ntimesteps 0 350]);
        %         hold off;
        %         subplot(3,2,2)
        %         title('Learned value');
        %         plot(tvec,v_func); xlim([-1 ntimesteps+1]);
        %         ylabel('expected value')
        % %         bar(c, mu_ij(i,:));
        % %         ylabel('basis function heights');
        %         %title('basis function values');
        %         %plot(tvec,v_jt);
        %         %ylabel('temporal basis function')
        % %         title(sprintf('trial # = %i', h)); %
        %
        %         subplot(3,2,3);
        %         scatter(rt_pred_i(1:i),rew_obs(1:i)); axis([1 ntimesteps 0 350]);
        % %         text(20, max(rew_obs), exptxt);
        %         hold on;
        %         plot(rt_pred_i(i),rew_obs(i),'r*','MarkerSize',20);  axis([1 ntimesteps 0 350]);
        %         hold off;
        %
        %         subplot(3,2,4);
        %         plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
        %         xlabel('time (centiseconds)')
        %         ylabel('uncertainty')
        %
        %         subplot(3,2,5);
        %         %eligibility trace
        %         title('eligibility trace');
        %         %elig_plot = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 1);
        %         %plot(tvec, elig_plot);
        %         plot(tvec, elig);
        %         xlabel('time(centiseconds)')
        %         ylabel('eligibility')
        %
        %
        %         subplot(3,2,6);
        %         plot(1:length(rt_obs), rt_obs, 'r');
        %         hold on;
        %         plot(1:length(rts_pred_exploit), rts_pred_exploit, 'b');
        % %        plot(1:length(rts_pred_explore), rts_pred_explore, 'k');
        % %        plot(1:length(rts_pred_explore), rts_uv_pred, 'g');
        %         hold off;
        %
        
        %% find unrewarded RTs
        unrew_rts = NaN(size(rt_obs));
        unrew_rts(rew_obs==0) = rt_obs(rew_obs==0);
        
        if strcmpi(modelname, 'value_softmax')
            %         subplot(2,1,1);
            %         plot(1:length(rt_obs), rt_obs, 'r');
            %         hold on;
            %         plot(1:length(rts_pred_exploit), rts_pred_exploit, 'b');
            %         hold off;
            %         subplot(2,1,2);
            %         contourf(1:ntrials, 1:ntimesteps, ret.v_it(1:ntrials,:)'); hold on;
            %         scatter(1:ntrials, rt_obs, 'r', 'Filled'); hold off;
            %         title('Value map');
            % %         pause(1);
        elseif strcmpi(modelname, 'uv')
            subplot(3,1,1);
            plot(1:length(rt_obs), rt_obs, 'r');
            hold on;
            plot(1:length(rts_pred_exploit), rts_pred_exploit, 'b');
            hold off;
            title('Red: actual RT, Blue: predicted RT');
            ax2 = subplot(3,1,2);
            contourf(1:ntrials, 1:ntimesteps, ret.v_it(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r','Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('Value map;   red: rewards,   blue: ommissions');
            colormap(ax2,summer);
            ax3 = subplot(3,1,3);
            contourf(1:ntrials, 1:ntimesteps, ret.uv(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r', 'Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('UV map; red: rewards,   blue: ommissions');
            colormap(ax3,summer);
        elseif strcmpi(modelname, 'v_discounted')
            subplot(3,1,1);
            plot(1:length(rt_obs), rt_obs, 'r');
            hold on;
            plot(1:length(rts_pred_exploit), rts_pred_exploit, 'b');
            plot(1:length(rts_pred_exploit_disc), rts_pred_exploit_disc, 'g*-');
            hold off;
            title('Red: actual RT, Blue: predicted RT_exploit, Green: predicted discounted RT_exploit');
            ax2 = subplot(3,1,2);
            contourf(1:ntrials, 1:ntimesteps, ret.v_it_undisc(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r','Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('Value map;   red: rewards,   blue: ommissions');
            colormap(ax2,summer);
            ax3 = subplot(3,1,3);
            contourf(1:ntrials, 1:ntimesteps, ret.v_it(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r', 'Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('Discounted value map; red: rewards,   blue: ommissions');
            colormap(ax3,summer);
            k = waitforbuttonpress;
        elseif strcmpi(modelname, 'uv_discounted')
            subplot(4,1,1);
            plot(1:length(rt_obs), rt_obs, 'r');
            hold on;
            plot(1:length(rts_pred_exploit), rts_pred_exploit, 'b');
            plot(1:length(rts_pred_exploit_disc), rts_pred_exploit_disc, 'g*-');
            plot(1:length(rts_pred_uv_disc), rts_pred_uv_disc, 'k');
            hold off;
            title('Red: actual RT, Blue: predicted RT_exploit, Green: predicted discounted RT_exploit, Black: discounted RT exploit with uncertainty premium');
            ax2 = subplot(4,1,2);
            contourf(1:ntrials, 1:ntimesteps, ret.v_it_undisc(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r','Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('Value map;   red: rewards,   blue: ommissions');
            colormap(ax2,summer);
            ax3 = subplot(4,1,3);
            contourf(1:ntrials, 1:ntimesteps, ret.v_it(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r', 'Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('Discounted value map; red: rewards,   blue: ommissions');
            colormap(ax3,summer);
            ax4 = subplot(4,1,4);
            contourf(1:ntrials, 1:ntimesteps, ret.uv(1:ntrials,:)'); hold on;
            scatter(1:ntrials, rt_obs,rew_obs+10, 'r', 'Filled');
            scatter(1:ntrials, unrew_rts,'b', 'Filled'); hold off;
            title('Discounted UV map; red: rewards,   blue: ommissions');
            colormap(ax4,summer);

            k = waitforbuttonpress;
            
            %         pause(3);
        end
        %figure(2); clf;
        %plot(tvec, u_func);
        %hold on;
        %plot(c, e_ij(i,:))
        %plot(c, e_ij(1:i,:)')
        %bar(c, sigma_ij(i,:))
        
        drawnow update;
        mov(i) = getframe(gcf);
    end
    
    %