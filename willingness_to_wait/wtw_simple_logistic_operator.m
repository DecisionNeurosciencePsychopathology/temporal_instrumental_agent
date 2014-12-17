%Author: Jonathan Wilson & Alex Dombrovski
%Date last modified: 10/8/14
%Matlab Version: R2012a
%Model of RT on the clock task using Gaussian stimuli for actions or "CSs"
%free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%learning rule - modified RW with TD(lambda)-type credit assignment
%choice rule: stochastic continuous RT weighted by value
%calls on RewFunction.m (Frank lab code)
%distrib_num = 1 - uniform, 2 - gen. Pareto, 3 - early beta, 4 - late beta


function [cost,constr,value_all,value_hist, wtw, mov] = wtw_simple_logistic_operator(params, distrib_num)
if (~exist('prr', 'var')), load prr; end;
%if (~exist('pev', 'var')), load pev; end;

number_of_stimuli = 12;
trial_plots = 1;


trial_length = 2000; %2000 10ms bins = 20s
min_trials = 100; % how many trials the agent will sample if he always waits
timeout = 80; % timeout after an immediate quit
task_length = min_trials*(trial_length+timeout);
maxtrials = ceil(task_length./timeout);          %maximum number of trials in a run
ntimesteps = trial_length; %number of timesteps in action space (set based on pseudorandom sampling of contingency above)
large_rew = 10; % the larger reward for waiting
small_rew = 1; % the smaller reward for quitting
sampling_reward = 0.01; % scaling constant representing the update to the
%uncertainty function in arbitrary units of uncertainty.  Scaled to make
%the sigmoid choice rule work.

%% initialize the move
mov=repmat(struct('cdata', [], 'colormap', []), maxtrials,1);

%%
pars = {{0, 1},{.4, .575, 0},{0.25, 0.25}, {5, 1}, {maxtrials,300}};
conds = {'unif', 'gp', 'beta', 'beta', 'normal'};
output_names = {'unif', 'gp', 'beta', 'late_beta', 'camel'};
%distrib_num = 4; % 1 - uniform, 2 - gen. Pareto, 3 - early beta, 4 - late
%beta, 5 - piecewise normal 

cond = output_names{distrib_num};
distrib_pars = pars{distrib_num};
distrib = conds{distrib_num};


constr = [];

if strcmpi(cond, 'unif')
    reward_times = prr.unif*ntimesteps;
%     ev_times = pev.unif*ntimesteps;
elseif strcmpi(cond, 'gp')
    reward_times = prr.gp*ntimesteps;
%     ev_times = pev.gp*ntimesteps;
elseif strcmpi(cond, 'beta')
    reward_times = prr.beta*ntimesteps;
%     ev_times = pev.beta*ntimesteps;
elseif strcmpi(cond, 'late_beta')
    reward_times = prr.late_beta*ntimesteps;
%     ev_times = pev.late_beta*ntimesteps;
elseif strcmpi(cond, 'camel')
    reward_times = prr.camel*ntimesteps;
    distrib = fitdist(prr.camel', 'kernel'); %Create PD object for EV computation
end


%--program start

%Initialize time step vector and allocate for memory
t=1:ntimesteps;


%Selected time points
step_num = number_of_stimuli - 1;
step = (ntimesteps-1)/step_num;
c = -step:step:ntimesteps+step;
nbasis=length(c);


%% free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
alpha = params(1);
lambda = params(2);
epsilon = params(3);
% log_k = params(4);
% k = exp(log_k);
% k is the hyperbolic discounting parameter
% sigma of Gaussian hard coded for now
sig = (150.*ntimesteps)./trial_length;


%% initialize RTs as a vector of zeros with a random initial RT;
wtw = zeros(1,maxtrials);

%Have initial rts be constant
%Willingness to wait (WTW) is the quit time, to which the model pre-commits
%at the start of the trial
%NB - times here are in 10ms (centisec) bins
%wtw(1) = ceil(.5*trial_length/10); <- trial_length/10 = deciseconds
wtw(1) = ceil(.5*trial_length);

% i = trial
% t = time step within trial, in centiseconds (1-500, representing 0-5 seconds)

% get a trial loop going
vh =            zeros(maxtrials, nbasis);     % initialize basis height values
deltah =        zeros(maxtrials, nbasis);     % prediction error assigned to each microstimulus
rh =            zeros(maxtrials, nbasis);     % reward assigned to each microstimulus
eh =            zeros(maxtrials, nbasis);     % eligibility traces for each microstimulus in relation to RT (US)
value_by_h =    zeros(nbasis, ntimesteps);   % value by microstimulus (rows for every microstimulus and columns for time points within trial)
value_hist =    zeros(maxtrials, ntimesteps); % history of value by trial
rew =           zeros(1, maxtrials);          % actual reward for each trial
%disc_rew =      zeros(1, maxtrials);          % discounted reward for each trial
uh =            zeros(maxtrials, nbasis);     % uncertainty at each microstimulus
udeltah =       zeros(maxtrials, nbasis);     % uncertainty prediction error for each microstumulus (hour)
sampling_h =    ones(maxtrials, nbasis);     % the amount of sampling assigned to each microstimulus at each trial
u_by_h =        zeros(nbasis, ntimesteps);   % uncertainty by microstimulus (rows for every microstimulus and columns for time points within trial)
u_hist =        zeros(maxtrials,ntimesteps);  % history of uncertainty by trial
ev     =        nan(1,maxtrials);             % expected rewards



reward_rate = zeros(1,maxtrials);
opportunity_cost = zeros(trial_length,maxtrials)';
cumulative_reward_fx = zeros(trial_length,maxtrials)';
return_on_policy = zeros(trial_length,maxtrials)';



%contruct gaussian matrix
gaussmat = zeros(nbasis,ntimesteps);

for j = 1:nbasis
    gaussmat(j,:) = gaussmf(t,[sig c(j)]);
end

temp_vh = zeros(nbasis,ntimesteps);
temp_uh = zeros(nbasis,ntimesteps);
gauss_params = [1.135 .99 .99 1.135];
corrective_basis = ones(1,ntimesteps);
corrective_basis = [gauss_params(1:2) ones(1,number_of_stimuli-2)...
gauss_params(3:4)]'*corrective_basis;


%Set up to run multiple runs for multiple ntrials
% for i = 1:ntrials
i=0;
task_time = 0;
while task_time<task_length;
    i=i+1;
    %         disp(i)
    if wtw(i)> reward_times(i)
        %% rew = real reward; disc_rew = reward discounted for the wait time
        rew(i) = large_rew;
        %disc_rew(i) = large_rew.*(1./(1+reward_times(i).*k));
        % get eligibility traces for each stimulus (h)
        % let's assume eligibility decays in inverse proporation to time
        % elapsed from or remaining until the peak of a stimulus
        %% NB - not sure if temporal generalization can distort the ...
        %contingency when the contingency is discontinuous, e.g. two neighboring bumps
        eh(i,:) = lambda.^(abs(reward_times(i)-c)+1);
        task_time = task_time+reward_times(i)+timeout;
        sampling_h(i,:) =  sampling_h(i,:)+sampling_reward; %sampling_h = vector of 0.01
    else
        rew(i) = small_rew;
        %disc_rew(i) = small_rew.*(1./(1+reward_times(i).*k));
        eh(i,:) = lambda.^(abs(wtw(i)-c)+1);
        task_time = task_time+wtw(i)+timeout;
        idx=find(eh(i,:)==max(eh(i,:)));
        sampling_h(i,idx:end) = (eh(i,idx:end)).*sampling_reward;
        sampling_h(i,1:idx) = max(eh(i,:)).*sampling_reward;
    end
    
    rh(i,:) = rew(i); 
    %rh(i,:) = disc_rew(i);
    
    %fprintf('i: %d, wtw(i): %.3f, reward_times(i): %.3f, task_time: %.3f\n', i, wtw(i), reward_times(i), task_time);
    
    %% stopped here on 10/9/14 at 6 PM.  Need cumulative sampling updates for wtw.
    %Sample cumulatively
    %idx=find(eh(i,:)==max(eh(i,:)));
    %sampling_h(i,idx:end) = (eh(i,idx:end));
    % learning rule: estimate value for each stimulus (h) separately
    deltah(i,:) = (rh(i,:) - vh(i,:)).*(eh(i,:));
    vh(i+1,:) = vh(i,:) + alpha.*deltah(i,:);
    
    %     udeltah(i,:) = uh(i,:) - sampling_h(i,:);
    %unlike value, uncertainty does not decay
    uh(i+1,:) = uh(i,:) - 0.1.*sampling_h(i,:); %0.1 is the fixed learning
    %rate for uncertainty updates
    
    %temp_vh = [gauss_params(1:2) ones(1,number_of_stimuli-2)...
    %gauss_params(3:4)]'*ones(1,ntimesteps);
    temp_vh = ones(1,ntimesteps);
    temp_vh = vh(i+1,:)'*temp_vh;
    temp_vh = corrective_basis.*temp_vh;
    
    value_by_h=temp_vh.*gaussmat(:,1:ntimesteps);
    
    
    % plot a sum of all stimulus value
    value_all = sum(value_by_h);
    
    
    temp_uh = ones(1,ntimesteps);
    temp_uh = uh(i+1,:)'*temp_uh;
    temp_uh = corrective_basis.*temp_uh;

    u_by_h=temp_uh.*gaussmat(:,1:ntimesteps);
    
    % plot a sum of uncertainty over all microstimuli
    u_all = sum(u_by_h);
    
    % store all the value functions by trial
    value_hist(i,:) = value_all;
    u_hist(i,:) = u_all;
    
    %% CHOICE RULE
    % find the RT corresponding to exploitative choice (choose randomly if value unknown)
    % NB: we added just a little bit of noise
    
    
    % 10/27/14 Updated choice Rule intergration/analytical method
    if i == 1   %indexing contingency
        reward_rate(i) = 0;
    else
        %update reward rate via delta rule
        reward_rate(i) = reward_rate(i-1) + alpha*(rew(i)/reward_times(i) - reward_rate(i-1));
    end
    
    %Opportunity cost = reward rate per trial * trial length
    opportunity_cost(i,:) = reward_rate(i).*(1:trial_length);
    cumulative_reward_fx(i,:) = cumtrapz(1:trial_length, value_all);
   
    
    return_on_policy(i,:) = cumulative_reward_fx(i,:)./opportunity_cost(i,:);
    
     if sum(value_all) == 0 || any(return_on_policy(i,:)==inf)
        %rt_exploit = ceil(rand(1)*ntrials); %random number within space
        rt_exploit = ceil(.5*trial_length); %default to mid-point of time domain
    else
        %rt_exploit = max(round(find(value_all==max(value_all))));
        rt_exploit = find(return_on_policy(i,:)==max(return_on_policy(i,:)));
        if rt_exploit > trial_length
            rt_exploit = trial_length;
        
        elseif rt_exploit < 0 %changed from if to elseif
            rt_exploit = 0;
        end
    end

  
%     if sum(value_all) == 0
%         %rt_exploit = ceil(rand(1)*ntrials); %random number within space
%         rt_exploit = ceil(.5*trial_length); %default to mid-point of time domain
%     else
%         %rt_exploit = max(round(find(value_all==max(value_all))));
%         rt_exploit = max(round(find(value_all(1:trial_length)==max(value_all))));%-5+round(rand(1).*10);
%         if rt_exploit > trial_length
%             rt_exploit = trial_length;
%         
%         elseif rt_exploit < 0 %changed from if to elseif
%             rt_exploit = 0;
%         end
%     end
    
    % find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)
    
    % u -- total amount of uncertainty on this trial (starts at 0 and decreases)
    u = mean(u_all);
    if u == 0
        rt_explore = ceil(.5*trial_length); %Changed!!! from 5000 previously
        
    else
        rt_explore = max(round(find(u_all(1:trial_length)==max(u_all))));
        
    end
    
    
    discrim = 100; %need to increase steepness of logistic given the tiny values we have here. Could free later
    sigmoid = 1/(1+exp(-discrim*(u - epsilon))); %Rasch model with epsilon as difficulty (location) parameter
    
    %hard classification of exploration for now at 0.5
    if i < maxtrials %do not populate rt on final trial
        if sigmoid > 0.5
            wtw(i+1) = rt_explore;
        else
            wtw(i+1) = rt_exploit;
        end
    end
    
    %% Compute the expected value of choice for the cost function
    if ~strcmpi(cond, 'camel')
        ev(i+1) = cdf(distrib, wtw(i+1)./trial_length, distrib_pars{:});
    else
        ev(i+1) = cdf(distrib, wtw(i+1)./trial_length);
    end
    
    if trial_plots == 1
        figure(1); %clf;
        
        %% make new intuitive plot like the one in Joe McGuire's papers
%         subplot(7,2,1:2)
%             title('black: wtw blue: RT(reward) red: RT(quit)');  hold on; ...
%             plot(find(rew(1:maxtrials)==large_rew),reward_times(rew(1:maxtrials)==large_rew),'bo-','LineWidth',2);
%         plot(find(rew(1:maxtrials)==small_rew),wtw(rew(1:maxtrials)==small_rew),'ro-','LineWidth',2); hold off;
%         
% 
%         subplot(7,2,3)
%         plot(t,value_all);
%         ylabel('value')
%         subplot(7,2,4)
%         %barh(sigmoid); axis([-.1 1.1 0 2]);
%         plot(t(1:ntimesteps),value_by_h);
%         ylabel('value temporal basis')
%         %title(sprintf('trial # = %i', h)); %
%         %         xlabel('time(ms)')
%         %         ylabel('reward value')
% 
%         subplot(7,2,5)
%         plot(cumulative_reward_fx(i,:));
%         ylabel('cumlulative reward function')
%         
%         subplot(7,2,6)
%         plot(return_on_policy(i,:));
%         ylabel('return on policy')
%         
%         %         subplot(7,2,5)
% %         plot(t(1:ntimesteps), u_all, 'r');
% %         xlabel('time (centiseconds)')
% %         ylabel('uncertainty')
% %         
% %         subplot(7,2,6)
% %         plot(t(1:ntimesteps),u_by_h);
% %         ylabel('uncertainty temporal basis')
% %         
%         subplot(7,2,7)
%         barh(sigmoid);
%         xlabel('explore or not'); axis([-.1 1.1 0 2]);
%         %         subplot(5,2,6)
%         %         barh(alpha), axis([0 .2 0 2]);
%         %         xlabel('learning rate')
%         subplot(7,2,8)
%         barh(lambda), axis([0.9 1 0 2]);
%         xlabel('decay')
%         
%         subplot(7,2,9)
%         barh(epsilon) %, axis([-0.5 0 0 2]);
%         xlabel('strategic exploration')
%         
%         subplot(7,2,10)
%         barh(u) %, axis([0 1000 0 2]);
%         xlabel('mean uncertainty')
%         %         pause(0.1);
%         
% %         subplot(7,2,11)
% %         plot(1:maxtrials,wtw)
% %         xlabel('will. to wait');
% %         
% %         subplot(7,2,12)
% %         barh(rt_exploit)
% %         xlabel('rt_exploit');
%         
%         
%         subplot(7,2,11) 
%         %barh(disc_rew(i));
%         barh(rew(i)), axis([0 10 1 1.099]);
%         ylabel('reward')
%         
%         subplot(7,2,12)
%         plot(1:i,find(value_all==max(value_all)));
%         xlabel('value_peak');
%         
%         
%         subplot(7,2,13:14) 
%         %plot(1:i,deltah(1:i));
%         barh(deltah(i,:), .5);
%         axis([-10 10 1 14])
%         xlabel('prediction error')
%        
%         
% %         xlabel('trial')
% %         ylabel('prediction error')
%         
%           % Plotexpected value for each trial
% %         subplot(6,2,12) 
% %         plot(1:i,ev(1:i));
% %         xlabel('trial')
% %         ylabel('ev')
%         
       
%% figures for the movie
subplot(4,2,1:4)    
title('black: wtw blue: RT(reward) red: RT(quit)');  hold on; ...
            plot(find(rew(1:maxtrials)==large_rew),reward_times(rew(1:maxtrials)==large_rew),'bo-','LineWidth',2);
        plot(find(rew(1:maxtrials)==small_rew),wtw(rew(1:maxtrials)==small_rew),'ro-','LineWidth',2); hold off;
        

        subplot(4,2,5)
        plot(t,value_all);
        ylabel('value')
        subplot(4,2,6)
        %barh(sigmoid); axis([-.1 1.1 0 2]);
        plot(t(1:ntimesteps),value_by_h);
        ylabel('value temporal basis')
        %title(sprintf('trial # = %i', h)); %
        %         xlabel('time(ms)')
        %         ylabel('reward value')

        subplot(4,2,7)
        plot(cumulative_reward_fx(i,:));
        ylabel('cumlulative reward function')
        
        subplot(4,2,8)
        plot(return_on_policy(i,:));
        ylabel('return on policy')

        drawnow update;
        mov(i) = getframe(gcf);
    end
    %     disp([i rts(i) rew(i) sum(value_all)])
end
%cost = -sum(rew);

cost = -sum(ev(~isnan(ev)));
%disp(cond);
%disp(params); disp(cost);
%     cond_matrix(j,:) = rts;
% end