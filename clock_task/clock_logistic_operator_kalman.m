%Kalman filter extension of temporal basis operator: each basis function tracks the mean value and uncertainty (sd) as a
%Gaussian.
%Authors: Jonathan Wilson, Alex Dombrovski, & Michael Hallquist
%Date last modified: 3/12/2015
%Matlab Version: R2012a
%Model of RT on the clock task using Gaussian stimuli for actions or "CSs"
%free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%learning rule - modified RW with TD(lambda)-type credit assignment
%choice rule: stochastic continuous RT weighted by value
%calls on RewFunction.m (Frank lab code)

function [cost,v_it,rts,mov,ret] = clock_logistic_operator_kalman(params, rngseeds, cond, ntrials, nbasis, ntimesteps, trial_plots, uvsum)
%params is the vector of parameters used to fit
%   params(1): alpha
%cond is the character string of the reward contingency
%ntrials is the number of trials to run
%nbasis is the number of radial basis functions used to estimate value and uncertainty
%ntimesteps is the number of time bins used for obtaining estimates of time functions for plotting etc.


%% free parameters: exploration tradeoff (epsilon), temporal decay (sd of Gaussian temporal spread)
epsilon = params(1); %NB, this is really c (u-v tradeoff) in uvsum model

%for now, setup the hack (just for prototype) that epsilon is the proportion reduction in variance at which
%exploration versus exploitation become equally probable in the sigmoid. Need to have a nicer choice rule that
%integrates time-varying information about value and uncertainty.

%note: Kalman filter does not have a free learning rate parameter.
if length(params) < 2
    disp('defaulting to prop_spread=.02');
    prop_spread = .02; %sigma of temporal spread function as a proportion of the finite interval (0..0.5)
else
    prop_spread = params(2);
end

if (length(params)) < 3
    disp('defaulting to no Gaussian random walk: k=0, s_grw=0');
    k=0;
    s_grw=0;
else
    k=params(3);
    s_grw=params(4);    
    %for now, don't link sigma for GRW to uncertainty (but could so that GRW decays with uncertainty)
end

if nargin < 2
    rngseeds=[98 83 66 10];
end
if nargin < 3, cond = 'DEV'; end
if nargin < 4, ntrials=100; end
if nargin < 5, nbasis = 24; end
if nargin < 6, ntimesteps=500; end
if nargin < 7, trial_plots = 1; end
if nargin < 8, uvsum = 0; end %whether to use uv sum for choice

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

%define radial basis
[c, sig, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(ntimesteps, nbasis, prop_spread);

%rescale s_grw wrt the interval (not as a proportion)
s_grw=s_grw*range(tvec); %determine SD of spread function

%add Gaussian noise with sigma = 1% of the range of the time interval to rt_explore
prop_expnoise=.01;
sig_expnoise=prop_expnoise*range(tvec);

%% initialize RTs chosen by agent as a nan vector;
rts = nan(1,ntrials);
%rts(1) = rand(1)*5000; %initial random rt?

%Have initial rts by random
%rts(1) = ceil(rand(1)*ntrials); %Changed!!! from 5000 previously

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
rew_i =         nan(1, ntrials);            % actual reward for each trial
ev_i =          nan(1, ntrials);            % expected value of choice for each trial (used for cost function)
u_jt =          zeros(nbasis, ntimesteps);  % uncertainty of each basis for each timestep
u_it =          zeros(ntrials,ntimesteps);  % history of uncertainty by trial at each timestep

%kalman setup
mu_ij =         nan(ntrials, nbasis);       %means of Gaussians for Kalman
k_ij =          zeros(ntrials, nbasis);     %Kalman gain (learning rate)
sigma_ij =      zeros(ntrials, nbasis);     %Standard deviations of Gaussians (uncertainty)

mu_ij(1,:) =    0; %expected reward on first trial is initialized to 0 for all Gaussians.

%noise in the reward signal: sigma_rew. In the Frank model, the squared SD of the Reward vector represents the noise in
%the reward signal, which is part of the Kalman gain. This provides a non-arbitrary initialization for sigma_rew, such
%that the variability in returns is known up front... This is implausible from an R-L perspective, but not a bad idea
%for getting a reasonable estimate of variability in returns. It would be interesting (plausible) to give a unique
%estimate to each basis based on its temporal receptive field, but hard to decide this up front. But it does seem like
%this could give local sensitivity to volatility (e.g., low- versus high-probability windows in time). For now, I'll
%just use the variance of the returns (ala Frank). But for the generative agent, this is not known up front -- so sample
%from the chosen contingency for each timestep as a guess.

%measurement noise
sigma_noise = repmat(std(arrayfun(@(x) RewFunction(x*10, cond, 0), tvec))^2, 1, nbasis);

%process noise should be some fixed proportion of measurement noise because
%k = U + 

u_threshold = (1-epsilon) * sigma_noise(1); %proportion reduction in variance from initial

%As in Frank, initialize estimate of std of each Gaussian to the noise of returns on a sample of the whole contingency.
%This leads to an effective learning rate of 0.5 since k = sigma_ij / sigma_ij + sigma_noise
sigma_ij(1,:) = sigma_noise; 

%fprintf('updating value by alpha: %.4f\n', alpha);
%fprintf('updating value by epsilon: %.4f with rngseeds: %s \n', epsilon, num2str(rngseeds));
fprintf('running agent with sigs: %.3f, epsilon: %.3f and rngseeds: %s \n', sig_spread, epsilon, num2str(rngseeds));

%10/16/2014 NB: Even though some matrices such as v_it have columns for discrete timesteps, the
%agent now learns entirely on a continuous time basis by maximizing the value and uncertainty curves over the
%radial basis set and estimating total uncertainty as the definite integral of the uncertainty function over the
%time window of the trial.

%rng('shuffle');

%objective expected value for this function
ev=[];
for val = 1:length(tvec)
    [~,ev(val)] = RewFunction(tvec(val).*10, cond);
end

figure(6); plot(tvec, ev);
title('Expected value of contingency');

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

    % estimate reward assigned to each stimulus h
    %[rew_i(i) ev_i(i)] = RewFunction(rts(i).*10, cond); %multiply by 10 because underlying functions range 0-5000ms
    %if (i < 100) 
    %    cond='IEV';
    %else
    %    cond='DEV';
    %end
    [rew_i(i) ev_i(i)] = RewFunction(rts(i).*10, cond); %multiply by 10 because underlying functions range 0-5000ms
            
    %Changing Kalman variance a posteriori should also use the elig*gain approach: [1 - k(ij)*elig(ij)]*sigma(ij)
    %this would only allow a 1.0 update*kalman gain for basis functions solidly in the window and a decay in diminishing
    %variance as the basis deviates from the timing of the obtained reward.

    %1) compute the Kalman gains for the current trial
    k_ij(i,:) = sigma_ij(i,:)./(sigma_ij(i,:) + sigma_noise);

    %2) update posterior variances on the basis of Kalman gains
    sigma_ij(i+1,:) = (1 - e_ij(i,:).*k_ij(i,:)).*sigma_ij(i,:);
    
    %3) update reward expectation
    delta_ij(i,:) = e_ij(i,:).*(rew_i(i) - mu_ij(i,:));
    mu_ij(i+1,:) = mu_ij(i,:) + k_ij(i,:).*delta_ij(i,:);
    
    v_jt=mu_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
    %v_jt=mu_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat_pdf; %use vector outer product to replicate weight vector
    
    %subjective value by timestep as a sum of all basis functions
    v_func = sum(v_jt);
    
    %uncertainty is now a function of Kalman uncertainties.
    u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat;
    %u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat_pdf;
    u_func = sum(u_jt);

    v_it(i+1,:) = v_func;
    u_it(i+1,:) = u_func;
    
    %% CHOICE RULE
    % find the RT corresponding to exploitative choice (choose randomly if value unknown)
    % NB: we added just a little bit of noise
    if sum(v_func) == 0
        %rt_exploit = ceil(rand(1)*ntrials); %random number within space
        rt_exploit = ceil(.5*ntimesteps); %default to mid-point of time domain
    else
        %DANGER: fminbnd is giving local minimum solution that clearly
        %misses the max value. Could switch to something more comprehensive
        %like rmsearch, but for now revert to discrete approximation
        %rt_exploit = fminbnd(@(x) -rbfeval(x, mu_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
        %figure(2);
        %vfunc = rbfeval(0:500, mu_ij(i+1,:), c, ones(1,nbasis).*sig);
        %plot(0:500, vfunc);
        rt_exploit = find(v_func==max(v_func));
        %rt_exploit = max(round(find(v_func==max(v_func(20:500)))))-5+round(rand(1).*10);
        if rt_exploit > max(tvec)
            rt_exploit = max(tvec);
        end
    end
        
    % find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)
    
    % u -- total amount of uncertainty on this trial (starts at 0 and decreases)
    % u = mean(u_func);
    
    %use integration to get area under curve of uncertainty
    %otherwise, our estimate of u is discretized, affecting cost function
    
    total_u = integral(@(x)rbfeval(x, sigma_ij(i+1,:), c, ones(1,nbasis).*sig), min(tvec), max(tvec));
    u = total_u/max(tvec); %make scaling similar to original sum? (come back)...
    
    if u == 0
        %rt_explore = ceil(rand(1)*ntrials);
        rt_explore = ceil(.5*ntimesteps); 
    else
        %rt_explore = fminbnd(@(x) -rbfeval(x, sigma_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
        rt_explore = find(u_func==max(u_func), 1) + round(sig_expnoise*randn(1,1)); %return position of first max and add gaussian noise
        %rt_explore = max(round(find(u_func==max(u_func(20:500)))));        
        
        %rt_explore=500; %always choose max for testing
        %rt_explore = randi([400,500],1);
    end

    %fprintf('trial: %d rt_exploit: %.2f rt_explore: %.2f\n', i, rt_exploit, rt_explore);
    
    discrim = 0.1;
    sigmoid = 1/(1+exp(-discrim.*(u - u_threshold))); %Rasch model with epsilon as difficulty (location) parameter
        
    %soft classification (explore in proportion to uncertainty)
    rng(explore_rng_state); %draw from explore/exploit rng
    choice_rand=rand;
    explore_rng_state=rng; %save state after random draw above
    
    %determine whether to strategic explore versus GRW
    rng(exptype_rng_seed);
    explore_type_rand=rand;
    exptype_rng_seed=rng;
    
    %determine step size for GRW
    rng(grw_step_rng_state); %draw from GRW rng
    grw_step=round(s_grw*randn(1,1));
    grw_step_rng_state=rng; %save state after draw 
    
    %rng('shuffle');
    if i < ntrials
        if choice_rand < sigmoid
            if (explore_type_rand > k)
                %strategic
                exptxt='strategic explore';
                rts(i+1) = rt_explore;
            else
                %grw
                exptxt='grw explore';
                rt_grw = rts(i) + grw_step;
                
                %N.B.: Need to have more reasonable GRW near the edge such that it doesn't just oversample min/max
                %e.g., perhaps reflect the GRW if rt(t-1) was already very close to edge and GRW samples in that
                %direction again.                
                if rt_grw > max(tvec), rt_grw = max(tvec);
                elseif rt_grw < min(tvec), rt_grw = min(tvec); end
                rts(i+1) = rt_grw;
            end
        else
            exptxt='exploit';%for graph annotation
            rts(i+1) = rt_exploit;
        end 
        
        %playing with basis update at the edge
        %rts(i+1) = randi([400,500],1); %force to late times
    
    end
    
    if uvsum == 1
        %new approach: use epsilon (actually c!) to scale relative contribution of u and v to bivariate choice
        %have left code above intact so that logistic model works as usual when uvsum==0
        uv=epsilon*v_func + (1-epsilon)*u_func;
        %uv=(epsilon*v_func) .* ((1-epsilon)*u_func);
        [~, rts(i+1)] = max(uv);
    end
    
    verbose=0;
    if verbose == 1
       fprintf('Trial: %d, Rew(i): %.2f, Rt(i): %.2f\n', i, rew_i(i), rts(i));
       fprintf('w_i,k:    '); fprintf('%.2f ', mu_ij(i,:)); fprintf('\n');
       fprintf('delta_ij:   '); fprintf('%.2f ', delta_ij(i,:)); fprintf('\n');
       fprintf('w_i+1,k:  '); fprintf('%.2f ', mu_ij(i+1,:)); fprintf('\n');
       fprintf('\n');
       
    end
    
    if trial_plots == 1
%         figure(1); clf;
%         subplot(5,2,1)
%         %plot(tvec,v_func);
%         scatter(rts(1:i),rew_i(1:i)); axis([1 500 0 350]);
%         hold on;
%         plot(rts(i),rew_i(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
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
%         plot(1:ntrials,rts, 'k');
%         ylabel('rt by trial'); axis([1 ntrials -5 505]);
        
        figure(1); clf;
        set(gca,'FontSize',18);
        subplot(3,2,1);
        title('Choice history');
        %plot(tvec,v_func);
        scatter(rts(1:i),rew_i(1:i)); axis([1 500 0 350]);
        text(20, max(rew_i), exptxt);
        hold on;
        plot(rts(i),rew_i(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
        hold off;
        subplot(3,2,2)
        title('Learned value');
        plot(tvec,v_func); xlim([-1 ntimesteps+1]);
        ylabel('expected value')
        subplot(3,2,3);
        
        %eligibility trace
        title('eligibility trace');
        %elig_plot = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 1);
        %plot(tvec, elig_plot);
        plot(tvec, elig);
%         bar(c, mu_ij(i,:));
%         ylabel('basis function heights');
        %title('basis function values');
        %plot(tvec,v_jt);
        %ylabel('temporal basis function')
%         title(sprintf('trial # = %i', h)); %
        xlabel('time(centiseconds)')
        ylabel('eligibility')
        
        subplot(3,2,4);
        plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
        xlabel('time (centiseconds)')
        ylabel('uncertainty')
        
        %figure(2); clf;
        %plot(tvec, u_func);
        %hold on;
        %plot(c, e_ij(i,:))
        %plot(c, e_ij(1:i,:)')
        %bar(c, sigma_ij(i,:))

        subplot(3,2,5);
        title('RT history');
        plot(1:ntrials, rts(1:ntrials));
        xlim([0 ntrials]);
        %plot(rts);
        ylabel('RT(ms)');
        
        
        if uvsum==1
            subplot(3,2,6);
            title('UV');
            plot(tvec,uv); xlim([-1 ntimesteps+1]);
            ylabel('UV')
        end

        
        drawnow update;
        mov(i) = getframe(gcf);
    end
    %     disp([i rts(i) rew_i(i) sum(v_func)])
end
%cost = -sum(rew_i);
cost = -sum(ev_i);

ret.mu_ij = mu_ij;
ret.sigma_ij = sigma_ij;
ret.k_ij = k_ij;
ret.delta_ij = delta_ij;
ret.e_ij = e_ij;
ret.v_it = v_it;
ret.u_it = u_it;
ret.rew_i = rew_i;
ret.ev_i = ev_i;


