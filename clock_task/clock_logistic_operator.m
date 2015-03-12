%Authors: Jonathan Wilson, Alex Dombrovski, & Michael Hallquist
%Date last modified: 3/12/2015
%Matlab Version: R2012a
%Model of RT on the clock task using Gaussian stimuli for actions or "CSs"
%free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%learning rule - modified RW with TD(lambda)-type credit assignment
%choice rule: stochastic continuous RT weighted by value
%calls on RewFunction.m (Frank lab code)

% function [cost,constr, value_func] = clock_smooth_action_model(...
%     params, cond, nbasis, trial_plots)
function [cost,constr,value_func,value_it,rts,mov] = clock_logistic_operator(params, rngseeds, cond, ntrials, nbasis, ntimesteps)
%params is the vector of parameters used to fit
%   params(1): alpha
%cond is the character string of the reward contingency
%ntrials is the number of trials to run
%nbasis is the number of radial basis functions used to estimate value and uncertainty
%ntimesteps is the number of time bins used for obtaining estimates of time functions for plotting etc.


%% free parameters: exploration tradeoff (epsilon), learning rate (alpha), temporal decay (sd of Gaussian temporal spread)
epsilon = params(1);

if length(params) < 2
    disp('defaulting to alpha=.06');
    alpha = .06;
else
    alpha = params(2);
end

if length(params) < 3
    disp('defaulting to prop_spread=.02');
    prop_spread = .02; %sigma of temporal spread function as a proportion of the finite interval (0..0.5)
else
    prop_spread = params(3);
end
    
if nargin < 2
    rngseeds=[98 83];
end
if nargin < 3, cond = 'DEV'; end
if nargin < 4, ntrials=100; end
if nargin < 5, nbasis = 24; end
if nargin < 6, ntimesteps=500; end

%states for random generators are shared across functions to allow for repeatable draws
global rew_rng_state explore_rng_state;

%initialize states for two repeatable random number generators using different seeds
rew_rng_seed=rngseeds(1);
explore_rng_seed=rngseeds(2);
rng(rew_rng_seed);
rew_rng_state=rng;
rng(explore_rng_seed);
explore_rng_state=rng;

trial_plots = 1; %whether to show trialwise graphics of parameters

%initialize movie storage
mov=repmat(struct('cdata', [], 'colormap', []), ntrials,1);

constr = [];

%Initialize time step vector and allocate for memory
tvec=1:ntimesteps;
sig_spread=prop_spread*range(tvec); %determine SD of spread function

%setup centers (means) and sds of basis functions
%based on testing in fix_rbf_basis.m, place the lowest center 12.5% below the first timestep and the last
%center 12.5% above last timestep. We were giving the agent 14 basis functions previously, so keep that for
%now (arg 6 above). SD should be calculated to give a Cohen's d of 1.52 between basis functions (~45% distribution overlap).

%margin_offset=0;
margin_offset = (max(tvec) - min(tvec))*.125; % 12.5% offset

%define lowest and highest centers
tmin = min(tvec) - margin_offset; tmax=max(tvec) + margin_offset;
c=tmin:(tmax-tmin)/(nbasis-1):tmax;

sig = (c(2) - c(1))/1.52; %cohen's d of 1.52 between basis functions

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
w_ij =          zeros(ntrials, nbasis);     % initialize basis height values
delta_ij =      zeros(ntrials, nbasis);     % prediction error assigned to each microstimulus
e_ij =          zeros(ntrials, nbasis);     % eligibility traces for each microstimulus in relation to RT (US)
value_jt =      zeros(nbasis, ntimesteps);  % value by microstimulus (rows for every microstimulus and columns for time points within trial)
value_it =      zeros(ntrials, ntimesteps); % history of value function by trial
rew_i =         nan(1, ntrials);            % actual reward for each trial
ev_i =          nan(1, ntrials);            % expected value of choice for each trial (used for cost function)
u_ij =          zeros(ntrials, nbasis);     % uncertainty at each microstimulus over trials (history)
u_jt =          zeros(nbasis, ntimesteps);  % uncertainty of each basis for each timestep
u_it =          zeros(ntrials,ntimesteps);  % history of uncertainty by trial at each timestep

%construct radial basis matrix using Gaussians
gaussmat = zeros(nbasis,ntimesteps);

for j = 1:nbasis
    gaussmat(j,:) = gaussmf(tvec,[sig c(j)]);
end

%normalize gauss functions to each have AUC = 1.0 within observed time interval
%this is essentially a truncated Gaussian basis such that AUC = 1.0 for basis functions within interval
maxauc=sum(gaussmat,2)*ones(1,length(tvec)); %outer product of vectors to allow for col-wise division below
gaussmat_trunc=gaussmat./maxauc;

%figure(20); plot(tvec,gaussmat); title('Regular RBF');
%figure(21); plot(tvec,gaussmat_trunc); title('Truncated RBF');

%NB: The working combination is: Regular Gaussian RBF for function evaluation, but weight update proceed by truncated
%Gaussian basis and truncated Gaussian spread function. Conceptually, this is like squeezing in all of the variance of
%the update inside the finite interval but then evaluating the function on the regular RBF... to be developed further
%since this remains a bit of a mystery. :)



%note: with truncated gaussian basis (which functions properly), the idea of extending the basis representation outside
%of the finite interval of interest does not apply (i.e., the functions are not defined outside of the bounds).
%gauss mat for all possible timesteps
%lowest_t = tmin - 4*sig; %3 SDs below lowest center.
%highest_t = tmax + 4*sig;

%t_all = round(lowest_t):round(highest_t);

%gaussmat_all = zeros(nbasis,length(t_all));

%for j = 1:nbasis
%    gaussmat_all(j,:) = gaussmf(t_all,[sig c(j)]);
%end
%plot(t_all, gaussmat_all);

%fprintf('updating value by alpha: %.4f\n', alpha);
%fprintf('updating value by epsilon: %.4f with rngseeds: %s \n', epsilon, num2str(rngseeds));
fprintf('running agent with alpha: %.3f, sigs: %.3f, epsilon: %.3f and rngseeds: %s \n', alpha, sig_spread, epsilon, num2str(rngseeds));

%10/16/2014 NB: Even though some matrices such as value_it have columns for discrete timesteps, the
%agent now learns entirely on a continuous time basis by maximizing the value and uncertainty curves over the
%radial basis set and estimating total uncertainty as the definite integral of the uncertainty function over the
%time window of the trial.

%rng('shuffle');

%Set up to run multiple runs for multiple ntrials
for i = 1:ntrials
    
    % get symmetric eligibility traces for each basis function (temporal generalization)
    % generate a truncated Gaussian basis function centered at the RT and with sigma equal to the free parameter.
    
    %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
    elig = gaussmf(tvec, [sig_spread, rts(i)]);
    
    %compute sum of area under the curve of the gaussian function
    auc=sum(elig);
    
    %divide gaussian update function by its sum so that AUC=1.0
    %note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
    %will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
    elig=elig/auc;
    
    %truncated gaussian eligibility
    %figure(7); plot(tvec, elig);
    
    %compute the intersection of the Gaussian spread function with the truncated Gaussian basis.
    %this is essentially summing the area under the curve of each truncated RBF weighted by the truncated
    %Gaussian spread function.
    e_ij(i,:) = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 2);

    % estimate reward assigned to each stimulus h
    [rew_i(i) ev_i(i)] = RewFunction(rts(i).*10, cond); %multiply by 10 because underlying functions range 0-5000ms
       
    %use the discrepancy between the expected value for each basis and the obtained reward
    %as prediction error and scale its effect on learning by eligibility.
    delta_ij(i,:) = e_ij(i,:).*(rew_i(i) - w_ij(i,:));
    w_ij(i+1,:) = w_ij(i,:) + alpha.*delta_ij(i,:);
    
    %lower lambdas will necessarily lead to a slower drop in average uncertainty because temporal precision 
    %of sampling is higher. this may be a problem for parameter fitting
    u_ij(i+1,:) = u_ij(i,:) - 1.0.*e_ij(i,:);
    
    value_jt=w_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
    
    %subjective value by timestep as a sum of all basis functions
    value_func = sum(value_jt);
    
    u_jt=u_ij(i+1,:)'*ones(1,ntimesteps) .*gaussmat;
    u_func = sum(u_jt);
        
    value_it(i,:) = value_func;
    u_it(i,:) = u_func;
    
    %% CHOICE RULE
    % find the RT corresponding to exploitative choice (choose randomly if value unknown)
    % NB: we added just a little bit of noise
    if sum(value_func) == 0
        %rt_exploit = ceil(rand(1)*ntrials); %random number within space
        rt_exploit = ceil(.5*ntimesteps); %default to mid-point of time domain
    else
        %DANGER: fminbnd is giving local minimum solution that clearly
        %misses the max value. Could switch to something more comprehensive
        %like rmsearch, but for now revert to discrete approximation
        %rt_exploit = fminbnd(@(x) -rbfeval(x, w_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
        %figure(2);
        %vfunc = rbfeval(0:500, w_ij(i+1,:), c, ones(1,nbasis).*sig);
        %plot(0:500, vfunc);
        rt_exploit = find(value_func==max(value_func));
        %rt_exploit = max(round(find(value_func==max(value_func(20:500)))))-5+round(rand(1).*10);
        if rt_exploit > 500
            rt_exploit = 500;
        end
    end
        
    % find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)
    
    % u -- total amount of uncertainty on this trial (starts at 0 and decreases)
    % u = mean(u_func);
    
    %use integration to get area under curve of uncertainty
    %otherwise, our estimate of u is discretized, affecting cost function
    
    total_u = integral(@(x)rbfeval(x, u_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
    u = total_u/500; %make scaling similar to original sum? (come back)...
    
    if u == 0
        %rt_explore = ceil(rand(1)*ntrials);
        rt_explore = ceil(.5*ntimesteps); 
    else
        %rt_explore = fminbnd(@(x) -rbfeval(x, u_ij(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
        rt_explore = find(u_func==max(u_func), 1); %return position of first max
        %rt_explore = max(round(find(u_func==max(u_func(20:500)))));        
        
        %rt_explore=500; %always choose max for testing
        %rt_explore = randi([400,500],1);
    end

    %fprintf('trial: %d rt_exploit: %.2f rt_explore: %.2f\n', i, rt_exploit, rt_explore);
    
    discrim = 50; %need to increase steepness of logistic given the tiny values we have here. Could free later
    sigmoid = 1/(1+exp(-discrim.*(u - epsilon))); %Rasch model with epsilon as difficulty (location) parameter
        
    %soft classification (explore in proportion to uncertainty)
    rng(explore_rng_state); %draw from reward rng
    %rng('shuffle');
    if i < ntrials
        rts(i+1) = rt_explore;
        if rand < sigmoid
            rts(i+1) = rt_explore;
        else
            rts(i+1) = rt_exploit;
        end 
        
        %playing with basis update at the edge
        %rts(i+1) = randi([400,500],1); %force to late times
    
    end
    explore_rng_state=rng; %save state after random draw above
    
    verbose=0;
    if verbose == 1
       fprintf('Trial: %d, Rew(i): %.2f, Rt(i): %.2f\n', i, rew_i(i), rts(i));
       fprintf('w_i,k:    '); fprintf('%.2f ', w_ij(i,:)); fprintf('\n');
       fprintf('delta_ij:   '); fprintf('%.2f ', alpha*delta_ij(i,:)); fprintf('\n');
       fprintf('w_i+1,k:  '); fprintf('%.2f ', w_ij(i+1,:)); fprintf('\n');
       fprintf('\n');
       
    end
    
    if trial_plots == 1
        figure(1); clf;
        subplot(5,2,1)
        %plot(tvec,value_func);
        scatter(rts(1:i),rew_i(1:i)); axis([1 500 0 350]);
        hold on;
        plot(rts(i),rew_i(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
        hold off;
        subplot(5,2,2)
        plot(tvec,value_func); xlim([-1 ntimesteps+1]);
        ylabel('value')
        subplot(5,2,3)
        
%         bar(c, w_ij(i,:));
%         ylabel('basis function heights');
        plot(tvec,value_jt);
        ylabel('temporal basis function')
%         title(sprintf('trial # = %i', h)); %
                xlabel('time(ms)')
                ylabel('reward value')
        
        subplot(5,2,4)
        plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
        xlabel('time (centiseconds)')
        ylabel('uncertainty')
        
        subplot(5,2,5)
        barh(sigmoid);
        xlabel('p(explore)'); axis([-.1 1.1 0 2]);
        subplot(5,2,6)
        barh(alpha), axis([0 .2 0 2]);
        xlabel('learning rate')
        subplot(5,2,7)
        barh(sig_spread), axis([0.0 0.5 0 2]);
        xlabel('decay')
        subplot(5,2,8)
        barh(epsilon), axis([-0.5 0 0 2]);
        xlabel('strategic exploration')
        
        subplot(5,2,9)
        barh(u) %, axis([0 1000 0 2]);
        xlabel('mean uncertainty')
        %         pause(0.1);
        subplot(5,2,10)
        plot(1:ntrials,rts, 'k');
        ylabel('rt by trial'); axis([1 ntrials -5 505]);
%         
%         figure(1); clf;
%         subplot(2,2,1)
%         %plot(tvec,value_func);
%         scatter(rts(1:i),rew_i(1:i)); axis([1 500 0 350]);
%         hold on;
%         plot(rts(i),rew_i(i),'r*','MarkerSize',20);  axis([1 500 0 350]);
%         hold off;
%         subplot(2,2,2)
%         plot(tvec,value_func); xlim([-1 ntimesteps+1]);
%         ylabel('value')
%         subplot(2,2,3)
%         
% %         bar(c, w_ij(i,:));
% %         ylabel('basis function heights');
%         plot(tvec,value_jt);
%         ylabel('temporal basis function')
% %         title(sprintf('trial # = %i', h)); %
%                 xlabel('time(ms)')
%                 ylabel('reward value')
%         
%         subplot(2,2,4)
%         plot(tvec, u_func, 'r'); xlim([-1 ntimesteps+1]);
%         xlabel('time (centiseconds)')
%         ylabel('uncertainty')
        
        %figure(2); clf;
        %plot(tvec, u_func);
        %hold on;
        %plot(c, e_ij(i,:))
        %plot(c, e_ij(1:i,:)')
        %bar(c, u_ij(i,:))

        drawnow update;
        mov(i) = getframe(gcf);
    end
    %     disp([i rts(i) rew_i(i) sum(value_func)])
end
%cost = -sum(rew_i);
cost = -sum(ev_i);

%     cond_matrix(j,:) = rts;
% end

% %% to be run after code completes
% for i = 1:length(cond_matrix)
%     plot(smooth(cond_matrix(i,:))) %Load in either iev or dev matrix and change variable name
%     pause(0.5)
% end



%figure(3)
% plot(smooth(t_max));
%clf;
%plot(mean(value_it));
% scatter(1:i,smooth(t_max),100);
% xlabel('trial')
% title(sprintf('time of maxumum value by trial'));
%figure(3)
%vh_tot = sum(vh);
%plot(1:4,vh_tot)
% tau_rt = rt;
% tau_max = c;
% tau_h = ones(5000,4); %Initalize
%
% %calculate deltaht, finding the distances between any rt in the time series
% %and the 4 cardnial time points, then dividing the singular rt point by the four distances.
% for i = 1:length(tau_rt)
%     tau_h(i,:)=abs(tau_rt(i) - tau_max);
%     deltaht(i,:) = rt(i)./tau_h(i,:);
% end

%% leftovers
% epsilon*u is the influence of uncertainty on choice
% since u tends to be negative, epsilon has to bring it to a sane range
%  epsilon/u.^2 is our weight;
%     u_weight = epsilon./(u^2+.001);
%         u=-.1;
%         epsilons = -100:1:100;
%         for ct = 1:length(epsilons)
%             epsilon = epsilons(ct);
%         sigmoid(ct) = 1./(1+exp((u.*epsilon)));
%         end
%         clf;
%         plot(epsilons,sigmoid);

% the choice probability (Moustafa 2008 after McClure 2003 after Egelman 1998)
% b and m are scaling constants
%         b = 2;
% need to scale m for p to add up to 1
%     m = 0.8;
%     p = (1./(1+exp(-m.*(value_func - b))));
%      if trial_plots == 1
%             subplot(3,1,3)
%
%         plot(tvec,p);
%         xlabel('time(ms)');
%         ylabel('choice probability');
%     pause(0.05)

% Now the actual choice rule
%     if sum(value_func) == 0
%         p_choice = (1./length(tvec-1)).*ones(1,ntimesteps);
%     else
%         p_choice = value_func./norm(value_func,1);
%
%     end
%     rts(i+1) = round(randsample(ntimesteps,1,true,p_choice));
%

%     for h = 1:nbasis
%         %plot stimuli h separately at trial level, has rows for every
%         %microstimulus and columns for time points within trial
%         value_jt(h,:) = vh(i,h)*gaussmf(tvec, [sig c(h)]);
%
%     end
    
% plot the location of the max(value)
%     if max(value_func) == 0
%         t_max(i) = (ntimesteps-1)./2;
%         else
%         t_max(i) = find(value_func == max(value_func));
%     end
% store all the value functions by trial

%get a vector of rewards for trial
%r = rand(1,100);
% r = zeros(1,ntrials);
% for i = 1:ntrials
%     r(i) = RewFunction(rts(i),cond);
% end

%Load in pseudo-random reward matrix
%pseudorand_rew_generator(ntrials,cond);

