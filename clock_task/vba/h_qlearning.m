function  [fx] = h_qlearning(x_t, theta, u, inF)
% This is the evolution function of the Q-learning model, it will output
% hidden state values for both waits and quits in a single column array.
% [fx,dfdx,dfdP] = h_qlearning(x,P,u,in)
%
% IN:
%   - x_t : Q values (100x1, [1:50] waits, [51:100] quits)
%   - theta : theta will be alpha and epsilon with fixed gamma and lambda
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: Q values (100x1, [1:50] waits, [51:100] quits)

%persistent Q;

%Set up free paramters -- epsilon might need to go in the ovservation function
alpha = 1./(1+exp(-theta(1)));
%epsilon = 1./(1+exp(-theta(2)));
gamma = .99;
lambda = .99;

%Set up variables from inputs
ntimesteps = inF.ntimesteps;
rt_obs = u(1); %essentially 1 = 100ms, 10 = 1000ms, 40 = 4000ms ect
rew = u(2);


% initialize Q with values from x_t
Q = [x_t(1:ntimesteps)  x_t(ntimesteps+1:end)];

%Set up eligibility traces
e = zeros(ntimesteps,ntimesteps);
for t = 1:ntimesteps
    for w = t:ntimesteps
        e(w,t)=(gamma^(w-t))*(lambda^(w-t));
    end
end


time_step = 1; %initialize current position to the start
next_time_step = 1; %initialize the next position to the start
% 1) choose initial action at current position using an epsilon greedy policy derived from Q
%a = action_list(1);

%Not needed right now
%[qmax, a] = chooseAction(epsilon, Q(time_step,:), time_step, ntimesteps);

episodeFinished = 0;
exit_criteria_1=0;
exit_criteria_2=0;
%Wander in the environment until agent 1) achieves the goal state, or 2) falls off the cliff.
%This constitutes a single episode (i.e., arrive at the final/absorbing state).

step = 1;
while(episodeFinished == 0 && time_step < ntimesteps)
    
    next_time_step = time_step+1; %next timestep
    
    %     if a == 2 %i.e., quit
    %         exit_criteria_2=1;
    %     end
    
    %Determine qmax for the next action, actions are predetermined based on
    %value from Q.
    %[qmax,a_next] = max(Q(next_time_step,:)); %qmax is the maximum value of Q, a is its index/position
    %[qmax, a_next] = chooseAction(epsilon, Q(next_time_step,:), next_time_step, ntimesteps);
    
    
    % update Q for prior location (t-1) based on payoff at this location (t)
    %curQ = Q(time_step, a); % working estimate of expected value at location(t-1), action(t-1)
    curQ = Q(time_step, 1); % working estimate of expected value at location(t-1), action(t-1)
    
    %Q-learning updates the curQ estimate based on the highest expected payoff for the next action (greedy)
    %This is an off policy decision because we have an epsilon greedy policy, not greedy
    %nextQ = qmax;
    nextQ = max(Q(next_time_step, :));
    
    %Update rule
    if time_step<= rt_obs
        %Quits
        if time_step==rt_obs
            %Added this since choose action is no longer used
            curQ = Q(time_step, 2); %Overwrite current value of curQ
            delta = rew-curQ;
            Q(time_step,2) = Q(time_step,2) + alpha*delta; %I had to code in the 1's and 2's for this part of the script
            if time_step==1
                %Q(:,a) = Q(:,a) + alpha*delta*e(time_step,:)'; %This was the line i was not sure about being correct
            else
                Q(:,1) = Q(:,1) + alpha*delta*e(time_step-1,:)';
            end
            exit_criteria_1=1;
            exit_criteria_2=1;
        %Waits
        else
            %delta = r+gamma*nextQ-curQ;
            delta = 0+gamma*nextQ-curQ;
            Q(:,1) = Q(:,1) + alpha*delta*e(time_step,:)';
            %Q(:,a) = Q(:,a) + alpha*delta*e(time_step,:)';
        end
    end
    
    %update time step
    time_step = next_time_step;

    %a = a_next;
    step = step + 1;
    
    
    %If the reward was harvested and the agent recorded a quit
    %we can finally exit the while loop
    if (exit_criteria_1 && exit_criteria_2)
        episodeFinished=1;
    end
end
fx = [Q(1:ntimesteps) Q(1+ntimesteps:end)]';

%fx = (fx-min(fx(:))) ./ (max(fx(:)-min(fx(:))));
% dfdx = fx * zeros(1,length(fx));
% dfdP = 0;

% function [qmax, a] = chooseAction(epsilon, Q, t, ntimesteps)
% %Function will return the max value at current timestep and the
% %appropriate action, based on value or a e-greedy policy
% %choose an action based on epsilon greedy strategy
% %does not permit selection of actions that move off the grid
%
% [qmax, a] =  max(Q); %qmax is the maximum value of Q, a is its index/position
% %
% % randomNormal = rand(1000,1);
% % if(randomNormal(randi(length(randomNormal),1)) <= epsilon)
% %     if rand>(1-(1/(ntimesteps-t)))
% %         a=2;
% %     else
% %         a=1;
% %     end
% % end
