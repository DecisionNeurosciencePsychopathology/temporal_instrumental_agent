%Authors: Jonathan Wilson, Alex Dombrovski, & Michael Hallquist
%Date last modified: 10/16/14
%Matlab Version: R2012a
%Model of RT on the clock task using Gaussian stimuli for actions or "CSs"
%free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%learning rule - modified RW with TD(lambda)-type credit assignment
%choice rule: stochastic continuous RT weighted by value
%calls on RewFunction.m (Frank lab code)

% function [cost,constr, value_all] = clock_smooth_action_model(...
%     params, cond, nbasis, trial_plots)
function [cost,constr,value_all,value_hist,rts] = clock_logistic_operator(params, rngseeds, cond, ntrials, nbasis, ntimesteps)
%params is the vector of parameters used in fitting
%cond is the character string of the reward contingency
%ntrials is the number of trials to run
%nbasis is the number of radial basis functions used to estimate value and uncertainty
%nevalpts is the number of time bins used for obtaining estimates of time functions for plotting etc.

if nargin < 2
    rngseeds=[98 83];
end
if nargin < 3, cond = 'IEV'; end
if nargin < 4, ntrials=100; end
if nargin < 5, nbasis = 12; end
if nargin < 6, ntimesteps=500; end

%old approach: load pseudorandom draws from each contingency
%we now use repeatable random numbers from rng in RewFunction.m
%if (~exist('b', 'var')), load b; end;
%if (~exist('b', 'var')), load pseudorand_trial; end;
%if (~exist('b', 'var')), load pseudorand_trial_5000; end;
% if strcmpi(cond, 'IEV')
%     r = b.iev_rew;
% elseif strcmpi(cond, 'DEV')
%     r = b.dev_rew;
% elseif strcmpi(cond, 'QUADUP')
%     r = b.quadup_rew;
% end
%ntimesteps = length(r); %number of timesteps in action space (set based on pseudorandom sampling of contingency above)
%ntimesteps = size(r, 2); %number of timesteps in action space (set based on pseudorandom sampling of contingency above)

global rew_rng_state explore_rng_state;

%generate states for two repeatable random number generators using different seeds
rew_rng_seed=rngseeds(1);
explore_rng_seed=rngseeds(2);
rng(rew_rng_seed);
rew_rng_state=rng;
rng(explore_rng_seed);
explore_rng_state=rng;

trial_plots = 0; %whether to show trialwise graphics of parameters

constr = [];

%scatter(b.rts, b.iev_rew)
%scatter(b.rts, b.dev_rew)

%--program start

%Initialize time step vector and allocate for memory
t=1:ntimesteps;

%Selected time points
step_num = nbasis - 1;
% step = (ntimesteps-1)/step_num;
% c = 0:step:ntimesteps-1;

%expand basis functions at edges
step = (ntimesteps-1)/step_num;
c = -step:step:ntimesteps+step;
nbasis=length(c);
% c = [0 1250 2500 3750 5000];

%Now set up the meat and potatoes (i.e. the learning rate code)
%The equation to set up is:
%Vht = gamma*Vht-1 + alpha*deltaht where:
%deltaht = rt/tau_h where:
%tau_h = abs(tau_t - tau_max)

%%NB: not sure how to use gamma (normally the discounting constant
%gamma = params(1); %Random number between 0 and 1
% gamma = 0.7;

%learning rate:
% alpha = .1; %params(2);
lambda = .95;
alpha=0.1;
%epsilon=-.06;

%% free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%alpha = params(1);
%lambda = params(2);
epsilon = params(1);
% sigma of Gaussian hard coded for now
sig = (350.*ntimesteps)./5000;


% % decay for values of stimuli h
% epsilon = params(3)

% rt = 1:5000;         %Any temporal point in the time series (excluding zero)

%get a vector of rts for each trial
% rts = rand(1,1000)*5000;

%% initialize RTs as a nan vector;
rts = nan(1,ntrials);
%rts(1) = rand(1)*5000; %initial random rt?

%Have initial rts by random
%rts(1) = ceil(rand(1)*ntrials); %Changed!!! from 5000 previously

%Have initial rts be constant for cost function
rts(1) = ceil(.5*ntimesteps);

% i = trial
% t = time step within trial, in centiseconds (1-500, representing 0-5 seconds)

vh =            zeros(ntrials, nbasis);     % initialize basis height values
deltah =        zeros(ntrials, nbasis);     % prediction error assigned to each microstimulus
rh =            zeros(ntrials, nbasis);     % reward assigned to each microstimulus
eh =            zeros(ntrials, nbasis);     % eligibility traces for each microstimulus in relation to RT (US)
value_by_h =    zeros(nbasis, ntimesteps);  % value by microstimulus (rows for every microstimulus and columns for time points within trial)
value_hist =    zeros(ntrials, ntimesteps); % history of value by trial
rew =           zeros(1, ntrials);          % actual reward for each trial
ev =            nan(1, ntrials);            % expected value of choice for each trial
uh =            zeros(ntrials, nbasis);     % uncertainty at each microstimulus
udeltah =       zeros(ntrials, nbasis);     % uncertainty prediction error for each microstumulus (hour)
sampling_h =    zeros(ntrials, nbasis);     % the amount of sampling assigned to each microstimulus at each trial
u_by_h =        zeros(nbasis, ntimesteps);  % uncertainty by microstimulus (rows for every microstimulus and columns for time points within trial)
u_hist =        zeros(ntrials,ntimesteps);  % history of uncertainty by trial

%construct gaussian matrix
gaussmat = zeros(nbasis,ntimesteps);

for j = 1:nbasis
    gaussmat(j,:) = gaussmf(t,[sig c(j)]);
end

plot(t,gaussmat);

% t_max = zeros(1,ntrials);

%fprintf('updating value by alpha: %.4f\n', alpha);
fprintf('updating value by epsilon: %.4f with rngseeds: %s \n', epsilon, num2str(rngseeds));

%QUESTION/CONCERN. Should value and uncertainty functions be summed over all plausible values, not just 0-500?
%This seems especially relevant for the integral of u... But then again, we are only using the AUC and the max
%value/uncertainty within the interval of interest, so maybe this doesn't matter whatsoever.

%10/16/2014 NB: Even though some matrices such as value_hist have columns for discrete timesteps, the
%agent now learns entirely on a continuous time basis by maximizing the value and uncertainty curves over the
%radial basis set and estimating total uncertainty as the definite integral of the uncertainty function over the
%time window of the trial.

%rng('shuffle');

%Set up to run multiple runs for multiple ntrials
for i = 1:ntrials
    %         disp(i)
    % get eligibility traces for each stimulus (h)
    % let's assume eligibility decays in inverse proporation to time
    % elapsed from or remaining until the peak of a stimulus
    eh(i,:) = lambda.^(abs(rts(i)-c)+1);
    %     eh(i,:) = 1./(abs(rts(i)-c)+1);
    
    % estimate reward assigned to each stimulus h
    %rew(i) = r(i, rts(i)); %RewFunction(rts(i).*10, cond); 
    [rew(i) ev(i)] = RewFunction(rts(i).*10, cond); 
    %[dummy ev(i)] = RewFunction(rts(i), cond); 
    
    %NB: This is problematic! Conceptually, we don't want to spread the point value of rewards in time and then
    %compare the spread reward to the expect reward at each time. This leads to big prediction errors far from
    %the chosen RT because only a fraction of the reward remains, which will almost always be worse than if temporally
    %distal reward were actually sampled.
    rh(i,:) = rew(i).*(eh(i,:));
    
    sampling_h(i,:) = (eh(i,:));
    % learning rule: estimate value for each stimulus (h) separately
    %deltah(i,:) = rh(i,:) - vh(i,:);
    
    %this seems more sensible: use the discrepancy between the expected value for each basis and the obtained reward
    %as predicition error and scale its effect on learning by eligibility.
    deltah(i,:) = eh(i,:).*(rew(i) - vh(i,:));
    vh(i+1,:) = vh(i,:) + alpha.*deltah(i,:);
    
    %udeltah(i,:) = uh(i,:) - sampling_h(i,:);
    %unlike value, uncertainty does not decay
    
    %lower lambdas will necessarily lead to a slower drop in average uncertainty because temporal precision of
    %sampling is higher.
    uh(i+1,:) = uh(i,:) - .01.*eh(i,:);
    
    temp_vh = vh(i+1,:)'*ones(1,ntimesteps);
    value_by_h=temp_vh.*gaussmat;
    
    %subjective value by timestep as a sum of all basis functions
    value_all = sum(value_by_h);
            
    temp_uh = uh(i+1,:)'*ones(1,ntimesteps);
    u_by_h=temp_uh.*gaussmat;
    u_all = sum(u_by_h);
    
    value_hist(i,:) = value_all;
    u_hist(i,:) = u_all;
    
    %% CHOICE RULE
    % find the RT corresponding to exploitative choice (choose randomly if value unknown)
    % NB: we added just a little bit of noise
    if sum(value_all) == 0
        %rt_exploit = ceil(rand(1)*ntrials); %random number within space
        rt_exploit = ceil(.5*ntimesteps); %default to mid-point of time domain
    else
        rt_exploit = fminbnd(@(x) -rbfeval(x, vh(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
        %rt_exploit = find(value_all==max(value_all));
        %rt_exploit = max(round(find(value_all==max(value_all(20:500)))))-5+round(rand(1).*10);
        if rt_exploit > 500
            rt_exploit = 500;
        end
    end
    
    % find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)
    
    % u -- total amount of uncertainty on this trial (starts at 0 and decreases)
    % u = mean(u_all);
    
    %use numerical integration to get area under curve of uncertainty
    %otherwise, our estimate of u is discretized, affecting cost function
    
    total_u = integral(@(x)rbfeval(x, uh(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
    u = total_u/500; %make scaling similar? (come back)...
    
    if u == 0
        %rt_explore = ceil(rand(1)*ntrials); %Changed!!! from 5000 previously
        rt_explore = ceil(.5*ntimesteps); %Changed!!! from 5000 previously        
    else
        rt_explore = fminbnd(@(x) -rbfeval(x, uh(i+1,:), c, ones(1,nbasis).*sig), 0, 500);
        %rt_explore = find(u_all==max(u_all));
        %rt_explore = max(round(find(u_all==max(u_all(20:500)))));        
    end

    discrim = 50; %need to increase steepness of logistic given the tiny values we have here. Could free later
    sigmoid = 1./(1+exp(-discrim.*(u - epsilon))); %Rasch model with epsilon as difficulty (location) parameter
    
    %hard classification of exploration for now at 0.5
%     if i < ntrials %do not populate rt on final trial
%         if sigmoid > 0.5
%             rts(i+1) = rt_explore;
%         else
%             rts(i+1) = rt_exploit;
%         end
%     end
    
    %soft classification (explore in proportion to uncertainty)
    rng(explore_rng_state); %draw from reward rng
    %rng('shuffle');
    if i < ntrials
        if rand < sigmoid %b.randdraws(i)
            rts(i+1) = rt_explore;
        else
            rts(i+1) = rt_exploit;
        end 
    end
    explore_rng_state=rng; %save state after random draw above
    
    if trial_plots == 1
        figure(1); clf;
        subplot(5,2,1)
        %plot(t,value_all);
        scatter(rts(1:i),rew(1:i)); axis([1 500 0 350]);
        subplot(5,2,2)
        plot(t,value_all); xlim([-1 ntimesteps+1]);
        ylabel('value')
        subplot(5,2,3)
        plot(t,value_by_h);
        ylabel('temporal basis function')
        %title(sprintf('trial # = %i', h)); %
        %         xlabel('time(ms)')
        %         ylabel('reward value')
        subplot(5,2,4)
        plot(t, u_all, 'r'); xlim([-1 ntimesteps+1]);
        xlabel('time (centiseconds)')
        ylabel('uncertainty')
        
        subplot(5,2,5)
        barh(sigmoid);
        xlabel('p(explore)'); axis([-.1 1.1 0 2]);
        subplot(5,2,6)
        barh(alpha), axis([0 .2 0 2]);
        xlabel('learning rate')
        subplot(5,2,7)
        barh(lambda), axis([0.9 1 0 2]);
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
        ylabel('rt by trial')
        
        drawnow update;
    end
    %     disp([i rts(i) rew(i) sum(value_all)])
end
%cost = -sum(rew);
cost = -sum(ev);

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
%plot(mean(value_hist));
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
%     p = (1./(1+exp(-m.*(value_all - b))));
%      if trial_plots == 1
%             subplot(3,1,3)
%
%         plot(t,p);
%         xlabel('time(ms)');
%         ylabel('choice probability');
%     pause(0.05)

% Now the actual choice rule
%     if sum(value_all) == 0
%         p_choice = (1./length(t-1)).*ones(1,ntimesteps);
%     else
%         p_choice = value_all./norm(value_all,1);
%
%     end
%     rts(i+1) = round(randsample(ntimesteps,1,true,p_choice));
%

%     for h = 1:nbasis
%         %plot stimuli h separately at trial level, has rows for every
%         %microstimulus and columns for time points within trial
%         value_by_h(h,:) = vh(i,h)*gaussmf(t, [sig c(h)]);
%
%     end
    
% plot the location of the max(value)
%     if max(value_all) == 0
%         t_max(i) = (ntimesteps-1)./2;
%         else
%         t_max(i) = find(value_all == max(value_all));
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

