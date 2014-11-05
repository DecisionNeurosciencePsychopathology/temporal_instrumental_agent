%Author: Jonathan Wilson & Alex Dombrovski
%Date last modified: 10/1/14
%Matlab Version: R2012a
%Model of RT on the clock task using Gaussian stimuli for actions or "CSs"
%free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%learning rule - modified RW with TD(lambda)-type credit assignment
%choice rule: stochastic continuous RT weighted by value
%calls on RewFunction.m (Frank lab code)

% function [cost,constr, value_all] = clock_smooth_action_model(...
%     params, cond, number_of_stimuli, trial_plots)
function [cost,constr,value_all,value_hist] = clock_discrete_action_model(params)

load b;

number_of_stimuli = 12;
trial_plots = 1;
cond = 'IEV';
%
% if nargin<2
%     number_of_stimuli = 12;
%     trial_plots = 0;
%     cond = 78;
% elseif nargin<3
%     number_of_stimuli = 12;
%     trial_plots = 0;
% end
constr = [];
if strcmpi(cond, 'IEV')
    r = b.iev_rew;
    %rts = pseudor_iev_rts;
elseif strcmpi(cond 'DEV')
    %r = fliplr(b.iev_rew);
    r = b.dev_rew;
end
trials = length(r);
%Initialize condition 12 = CEV; 34 = CEVR; 56 = DEV; 78 = IEV
% cond = 56;



%scatter(b.rts, b.iev_rew)
%scatter(b.rts, b.dev_rew)



% trial_plots = 1;

%--program start

%Initialize time scale and allocate for memory
%t=0:5000;
%t=0:trials;
t=1:5000; %<-changed this on 10/2/14
%y=zeros(4,length(t));
%value_by_h = zeros(4,5001);
%value_all = zeros(1,5001);

%Selected time points
% c = [0 1250 2500 3750 5000];
step_num = number_of_stimuli;
step = (length(t)-1)/step_num;
c = 0:step:length(t)-1;


% c = [0 1250 2500 3750 5000];


% for i = 1:size(y,1)
%     y(i,:)=2*gaussmf(t,[sig c(i)]);
% end
% plot(t,y)
% xlabel('time(ms)')

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
% lambda = .99;

%% free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
alpha = params(1);
%lambda = params(2);
epsilon = params(2);
% sig hard coded for now
%sig = 500;
%sig = (250.*length(r))./5000;


% % decay for values of stimuli h
% epsilon = params(3)

% rt = 1:5000;         %Any temporal point in the time series (excluding zero)

%get a vector of rts for each trial
% rts = rand(1,1000)*5000;

%% initialize RTs as a vector of zeros with a random initial RT;
rts = zeros(1,trials);
%rts(1) = rand(1)*5000;

%Have initial rts by random
%rts(1) = ceil(rand(1)*trials); %Changed!!! from 5000 previously

%Have initial rts be constant
rts(1) = ceil(.5*trials);

%plot (1:100,rts)

%get a vector of rewards for trial
%r = rand(1,100);
% r = zeros(1,length(rts));
% for i = 1:length(rts)
%     r(i) = RewFunction(rts(i),cond);
% end

%Load in pseudo-random reward matrix
%pseudorand_rew_generator(trials,cond);

% t = time step, in ms
% i = trial

% get a trial loop going
%initialize values
vh(1,:) = zeros(1,length(c));

% deltah -- prediction error assigned to each microstimulus
deltah = zeros(length(rts),length(c));

% rh -- reward assigned to each microstimulus
rh = zeros(length(rts),length(c));

% eh - eligibility traces for each microstimulus in relation to RT (US)
eh = zeros(length(rts),length(c));

%value by microstimulus (rows for every microstimulus and columns for time
%points within trial)
value_by_h = zeros(length(c),length(t));

% value_hist - history of value by trial
value_hist = zeros(length(rts),length(t));

% rew -- actual reward for the trial
rew = zeros(1, length(r));

%%
%uh -- uncertainty at each microstimulus
uh = zeros(length(rts),length(c));

%u_deltah -- uncertainty prediction error for each microstumulus (hour)
%udeltah = zeros(length(rts),length(c));

% sampling_h -- the amount of sampling assigned to each microstimulus at
% each trial
sampling_h = zeros(length(rts),length(c));

%uncertainty by microstimulus (rows for every microstimulus and columns for time
%points within trial)
u_by_h = zeros(length(c),length(t));

% u_hist - history of uncertainty by trial
u_hist = zeros(length(rts),length(t));


% t_max = zeros(1,length(rts));
% runs = 100;

%Set up to run multiple runs for multiple trials
% for j = 1:runs
for i = 1:length(rts) %-1 %-1 % added -1 here!!!!!
    %         disp(i)
    % get eligibility traces for each stimulus (h)
    % let's assume eligibility decays in inverse proporation to time
    % elapsed from or remaining until the peak of a stimulus
    eh(i,:) = rts(i)
    
    
    rtdiffs=rts(i) - c;
    guyFed=rtdiffs == min(rtdiffs(find(rtdiffs > 0))); %left edge of stimulus receiving update
    eh(i,guyFed) = 1;

    %% stopped here 4pm 7Oct2014
    rew(i) = r(rts(i));
    rh(i,:) = rew(i).*(eh(i,:));
    
    sampling_h(i,:) = (eh(i,:));
    % learning rule: estimate value for each stimulus (h) separately
    deltah(i,:) = rh(i,:) - vh(i,:);
    vh(i+1,:) = vh(i,:) + alpha.*deltah(i,:);
    
    
%     udeltah(i,:) = uh(i,:) - sampling_h(i,:);
    %unlike value, uncertainty does not decay
    uh(i+1,:) = uh(i,:) - alpha.*eh(i,:);
    
    %clf
    
    for h = 1:length(c)
        %plot stimuli h separately at trial level, has rows for every
        %microstimulus and columns for time points within trial
        value_by_h(h,:) = vh(i,h)*gaussmf(t, [sig c(h)]);
        
    end
    
    % plot a sum of all stimulus value
    value_all = sum(value_by_h);
    
    for h = 1:length(c)
        %plot stimuli h separately
        u_by_h(h,:) = uh(i,h)*gaussmf(t, [sig c(h)]);
        
    end
    
    % plot a sum of uncertainty over all microstimuli
    u_all = sum(u_by_h);
    
    
    %         pause(0.05)
    
    
    % plot the location of the max(value)
    %     if max(value_all) == 0
    %         t_max(i) = (length(t)-1)./2;
    %         else
    %         t_max(i) = find(value_all == max(value_all));
    %     end
    % store all the value functions by trial
    value_hist(i,:) = value_all;
    u_hist(i,:) = u_all;
    
    %% the choice probability (Moustafa 2008 after McClure 2003 after Egelman 1998)
    % b and m are scaling constants
    %         b = 2;
    %% need to scale m for p to add up to 1
    %     m = 0.8;
    %     p = (1./(1+exp(-m.*(value_all - b))));
    %      if trial_plots == 1
    %             subplot(3,1,3)
    %
    %         plot(t,p);
    %         xlabel('time(ms)');
    %         ylabel('choice probability');
    %     pause(0.05)
    %%
    % Now the actual choice rule
    %     if sum(value_all) == 0
    %         p_choice = (1./length(t-1)).*ones(1,length(t));
    %     else
    %         p_choice = value_all./norm(value_all,1);
    %
    %     end
    %     rts(i+1) = round(randsample(length(t),1,true,p_choice));
    %
    
    %% CHOICE RULE
    % find the RT corresponding to exploitative choice (choose randomly if
    % value unknown)
    if sum(value_all) == 0
%         rt_exploit = ceil(rand(1)*trials); %Changed!!! from 5000 previously
        rt_exploit = ceil(.5*trials); %Changed!!! from 5000 previously

    else
        rt_exploit = max(round(find(value_all==max(value_all))));
    end
    
    % find the RT corresponding to uncertainty-driven exploration (try
    % random exploration if uncertainty is uniform)
    
    % u -- total amount of uncertainty on this trial (starts at 0 and
    % decreases)
    u = mean(u_all);
    if u == 0
%         rt_explore = ceil(rand(1)*trials); %Changed!!! from 5000 previously
        rt_explore = ceil(.5*trials); %Changed!!! from 5000 previously

    else
        rt_explore = max(round(find(u_all==max(u_all))));
    end
    
    % epsilon*u is the influence of uncertainty on choice
    % since u tends to be negative, epsilon has to bring it to a sane range
    %  epsilon/u.^2 is our weight; 
    u_weight = epsilon./(u^2+.001);
    if i<trials
    rts(i+1)= round(wmean([rt_exploit rt_explore],[1 u_weight],2)); 
    end
        
    %% also try logistic between max(value) and max(uncertainty)
    
    if trial_plots == 1
        figure(1); clf;
        subplot(5,2,1)
        %plot(t,value_all);
        scatter(rts(1:trials),rew);axis([1 500 1 350]);
        subplot(5,2,2)
        plot(t,value_all);
        ylabel('value')
        subplot(5,2,3)
        plot(t,value_by_h);
        ylabel('temporal basis function')
        %title(sprintf('trial # = %i', h)); %
%         xlabel('time(ms)')
%         ylabel('reward value')
        subplot(5,2,4)
        plot(t, u_all, 'r');
        xlabel('time(ms)')
        ylabel('uncertainty')
        subplot(5,2,5)
        barh(u_weight,'r')
        xlabel('uncertainty weight')
        subplot(5,2,6)
        barh(alpha), axis([0 .5 0 2]);
        xlabel('learning rate')
        subplot(5,2,7)
        barh(lambda), axis([0.9 1 0 2]);
        xlabel('decay')
        subplot(5,2,8)
        barh(epsilon), axis([0 100 0 2]);
        xlabel('strategic exploration')
        
        subplot(5,2,9)
        barh(u) %, axis([0 1000 0 2]);
        xlabel('mean uncertainty')
%         pause(0.1);
        subplot(5,2,10)
        plot(1:500,rts, 'k');
        ylabel('rt by trial') 
    end
%     disp([i rts(i) rew(i) sum(value_all)])
end
    cost = -sum(rew);

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

%--end of program