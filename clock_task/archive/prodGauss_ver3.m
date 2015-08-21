%Author: Jonathan Wilson
%Date last modified: 9/29/14
%Matlab Version: R2012a
%This is a script to generate gaussian distributions along a time axis of
%5000ms....so far.

trial_plots = 1;
%--program start

%Initialize time scale and allocate for memory
%t=0:5000;
t = 0:500;
y=zeros(4,length(t));
%value_by_h = zeros(4,5001);
%value_all = zeros(1,5001);

%Selected time points
% c = [0 1250 2500 3750 5000];
step_num = 10;
step = (length(t)-1)/step_num;
c = 0:step:length(t)-1;
% c = [0 1250 2500 3750 5000];


% sig hard coded for now
sig = 50;

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


%learning rate:
alpha = .1; %params(2);
lambda = .99;

% % decay for values of stimuli h
% epsilon = params(3)

% rt = 1:5000;         %Any temporal point in the time series (excluding zero)

%Initialize condition 12 = CEV; 34 = CEVR; 56 = DEV; 78 = IEV
cond = 'IEV';

%get a vector of rts for each trial
% rts = rand(1,1000)*5000;

%% initialize RTs as a vector of zeros with a random initial RT;
rts = zeros(1,100);
rts(1) = rand(1)*5000;
%plot (1:100,rts)

%get a vector of rewards for trial
%r = rand(1,100);
r = zeros(1,length(rts));
for i = 1:length(rts)
    r(i) = RewFunction(rts(i),cond);
end
% t = time step, in ms
% i = trial

% get a trial loop going
%initialize values
vh(1,:) = zeros(1,length(c));
deltah = zeros(length(rts),length(c));
rh = zeros(length(rts),length(c));
eh = zeros(length(rts),length(c));
value_by_h = zeros(length(c),length(t));
value_hist = zeros(length(rts),length(t));
t_max = zeros(1,length(rts));
runs = 100;

%Set up to run multiple runs for multiple trials
% for j = 1:runs
    for i = 1:length(rts) %-1 % added -1 here!!!!!
        disp(i)
        % get eligibility traces for each stimulus (h)
        % let's assume eligibility decays in inverse proporation to time
        % elapsed from or remaining until the peak of a stimulus
        eh(i,:) = lambda.^(abs(rts(i)-c)+1);
        %     eh(i,:) = 1./(abs(rts(i)-c)+1);
        
        % estimate reward assigned to each stimulus h
        %rh(i,:) = r(i).*(eh(i,:));
        rh(i,:) = r(i).*(eh(i,:));
        
        % learning rule: estimate value for each stimulus (h) separately
        deltah(i,:) = rh(i,:) - vh(i,:);
        vh(i+1,:) = (1-alpha).*vh(i,:) + alpha.*deltah(i,:);
        %vh(i+1,:) = vh(i,:) + alpha.*deltah(i,:);
        
        clf
        
        for h = 1:length(c)
            %plot stimuli h separately
            value_by_h(h,:) = vh(i,h)*gaussmf(t, [sig c(h)]);
            
            if trial_plots == 1
                figure(1); %clf;
                subplot(3,1,1)
                plot(t,value_by_h);
                %title(sprintf('trial # = %i', h)); %
                xlabel('time(ms)')
                ylabel('reward value')
            else
            end
            
        end
        % plot a sum of all stimulus value
        value_all = sum(value_by_h);
        if trial_plots == 1
            
            subplot(3,1,2)
            
            %plot(t,value_all);
            
            plotyy(t,value_all,t,r(i));
            pause(0.05)
        else
        end
        % plot the location of the max(value)
        %     if max(value_all) == 0
        %         t_max(i) = (length(t)-1)./2;
        %         else
        %         t_max(i) = find(value_all == max(value_all));
        %     end
        % store all the value functions by trial
        value_hist(i,:) = value_all;
        
        %% the choice probability (Moustafa 2008 after McClure 2003 after Egelman 1998)
        % b and m are scaling constants
        b = 2;
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
        if sum(value_all) == 0
            p_choice = (1./length(t-1)).*ones(1,length(t));
        else
            p_choice = value_all./norm(value_all,1);
        end
        rts(i+1) = randsample(length(t),1,true,p_choice);
        
        
    end
    
%     cond_matrix(j,:) = rts;
% end

%% to be run after code completes
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