function [maxReward, constr, quits,cumReward] = ClockWalking_3D(options)
global episodesStruct;
% CliffWalking: implements the cliff-walking problem using either the Q-Learning method or SARSA.

% You can pass the parameters of the problem using the options
% structure. Otherwise, default settings are used for running the program.

% written by: Sina Iravanian - June 2009
% sina@sinairv.com
% Please send your comments or bug reports to the above email address.

% adapted by: Michael Hallquist - February 2013

% Notes: The SARSA name reflects the fact that the main function for updating the Q-value depends on the 
%    current state of the agent "S1", the action the agent chooses "A1", the reward "R" the agent gets for
%    choosing this action, the state "S2" that the agent will now be in after taking that action, and finally
%    the next action "A2" the agent will choose in its new state.
%    Taking every letter in the quintuple (st, at, rt, st+1, at+1) yields the word SARSA.
%
%    Q-Learning versus SARSA: The difference may be explained as SARSA learns the Q values associated with 
%    taking the policy it follows itself, while Watkin's Q-learning learns the Q values associated with taking
%    the exploitation policy while following an exploration/exploitation policy.
%    Source: http://en.wikipedia.org/wiki/SARSA

% simps won't run without constr
constr=[];

if(nargin < 1),
    gamma = 0.9; %discount factor: how important are are short-term versus long-term gains (0 = short-sighted, 1 = long-term payoff)
    alpha = 0.05; %learning rate: how much to weight new information (0 = no learning, 1 = only value latest information)
    epsilon = 0.1; %how often does agent explore? 10% here
    ntimesteps = 50; %how many timesteps for the trial
    fontsize = 5;
    showTitle = 1;

    episodeCount = 2000;
    plotEpisodes = [20 200 700 1000 2000];
    allowDiagonal = 0; 
    canHold = 0;

    start = 1; %always start at first timestep
    agent = 'qlearning';
    speedbumpcol = 0; %by default, leave out speedbump
    speedbumpcost = -10;
    episodesToTrack = [ 5 100 200 300 1000 2000 ]; %which episodes to track (snapshots)
elseif ~isstruct(options) %This is only for running simps
    load('clock_options.mat')
    alpha=options(1);
    gamma=options(2);
    epsilon=options(3);
    
    ntimesteps = clock_options.ntimesteps;
    agent = clock_options.agent;
    episodeCount = 1000; %Split 500 IEV 500 DEV
    start = clock_options.start;
    episodesToTrack = clock_options.snapshotEpisodes;
    %cond = clock_options.cond;
    diagnos = 0;
    decay_flag = 0;
else
    gamma = options.gamma;
    alpha = options.alpha;
    epsilon = options.epsilon;
    ntimesteps = options.ntimesteps;
    fontsize = options.fontsize;
    showTitle = options.showTitle;

    episodeCount = options.episodeCount;
    plotEpisodes = options.plotEpisodes;
    allowDiagonal = options.allowDiagonal; 
    canHold = options.canHold;

    start = options.start;
    agent = options.agent;
    speedbumpcol = options.speedbumpcol;
    speedbumpcost = options.speedbumpcost;
    episodesToTrack = options.snapshotEpisodes;
    smallreward = options.smallreward;
    
    cond = options.cond;
    diagnos = options.diag;
    decay_flag = options.decayflag;
    wait_punishment = options.waitflag;
end

%states for random generators are shared across functions to allow for repeatable draws
global rew_rng_state explore_rng_state;
rngseeds=[98 83];
%initialize states for two repeatable random number generators using different seeds
rew_rng_seed=rngseeds(1);
explore_rng_seed=rngseeds(2);
rng(rew_rng_seed);
rew_rng_state=rng;
rng(explore_rng_seed);
explore_rng_state=rng;

%The amount of decay you want to decrease epsilon by
decay_val=.95; %Currently 5%



selectedEpIndex = 1;
%quit or wait
nactions = 2;

% initialize Q with zeros (no initial information about environment)
% Q is the action-value function: expected reward for each possible action at a given timestep in the trial.
global Q;
Q = zeros(ntimesteps, nactions, nactions);
cumReward = zeros(episodeCount, 1); %vector to track total reward earned per episode

a = 0; % an invalid action
% loop through episodes
% agent is allowed to sample environment many times and carry forward its learning to subsequent episodes.

%To show the updates in Q as the agent moves through environment, take full snapshots of some episodes
e2t_i = 1; %index of episode being tracked
episodesStruct = []; %will be a vector (episodes) of structs containing cell arrays of Q and a (per action/step)

%set factor depending upon a 10ms or 100ms bin size
factor=100; %bins are 100ms

quits=zeros(episodeCount,1);
explore_count=zeros(ntimesteps,1);


for ei = 1:episodeCount,
    
    if ~isstruct(options)
        if ei<=500 %Run 500 trials of IEV then 500 trials of DEV
            cond='IEV';
        else
            cond='DEV';
        end
    end
    
    %commented out for simps
    %fprintf('Running episode %d\n', ei);
    
    curpos = start; %initialize current position to the start
    nextpos = start; %initialize the next position to the start
    
    %decay epsilon
    if (decay_flag ==1 && episodeCount ~=1),epsilon=epsilon.*decay_val; end
    
    % 1) choose initial action at current position using an epsilon greedy policy derived from Q
    [qmax, a, explore] = chooseAction(epsilon, Q, curpos, ntimesteps, nactions);
      
    episodeFinished = 0;
    %Wander in the environment until agent 1) achieves the goal state, or 2) falls off the cliff.
    %This constitutes a single episode (i.e., arrive at the final/absorbing state).
        
    step = 1;
    while(episodeFinished == 0 && curpos.row < ntimesteps)

        
        
        % 2) take action a to move from curpos to nextpos
        % 3) receive reward r (at nextpos)
        nextpos.row = curpos.row+1; %next timestep
        %fprintf('ei is: %d, nextpos is %d\n', ei, nextpos);
        %fprintf('Location: %d %d\n', nextpos, nextpos);
        if a == 2 %i.e., quit
            %if we have arrived at the goal state, mark as finished and provide reward 0
            episodeFinished = 1;
            r = RewFunction(nextpos.row.*factor, cond);
            %fprintf('harvested rew: %.2f ', r);
            quits(ei)=curpos.row;
            
            
        else %wait
            r = 0;
%             curpos = curpos + 1;
%             [qmax, a, explore] = chooseAction(epsilon, Q, curpos, ntimesteps, nactions);
%             continue;
            %r = RewFunction(nextpos.*factor, cond); %every move recives probablistic reward
        end
        
        %if requested, take snapshot of agent in episode before next action chosen and Q is updated.
        if (e2t_i <= length(episodesToTrack) && ei == episodesToTrack(e2t_i))
            %Use cell array within each struct to allow for variable numbers of trials
            episodesStruct(e2t_i).Q{step} = Q; %state-value function (expected value)
            episodesStruct(e2t_i).a{step} = a; %action chosen
            episodesStruct(e2t_i).r{step} = r; %reward experienced
            episodesStruct(e2t_i).explore{step} = explore;
        end
        
        % 3) choose next action (a_next) based on maximum expected payoff at the new position: Q(nextpos)
        % use epsilon greedy policy based on Q
        
        % exploit: choose best action based on max Q (greedy selection)
        % put this here (not within explore/exploit if statement) so that qmax is always known for Q-learning.
        [qmax, a_next, explore] = chooseAction(epsilon, Q, nextpos, ntimesteps, nactions);
        
        % update Q for prior location (t-1) based on payoff at this location (t)
        curQ = Q(curpos.row, curpos.col, a); % working estimate of expected value at location(t-1), action(t-1) 
        if strcmp(agent, 'qlearning')
            %Q-learning updates the curQ estimate based on the highest expected payoff for the next action (greedy)
            %This is an off policy decision because we have an epsilon greedy policy, not greedy
            nextQ = qmax;
        elseif strcmp(agent, 'sarsa')
            %SARSA updates curQ using information about expected payoff at the next time step, knowing the next action
            %This is on policy in this sense that it respects whatever action selection policy is used
            %  (here, epsilon-greedy instead of pure greedy)
            nextQ = Q(nextpos.row, nextpos.col, a_next);
        end

        % update Q matrix (action-value)
        % learned value is [r + gamma*nextQ]
        % old value is curQ
        
        %I think this is what Alex wanted, once the agent has decided it
        %wants to go into the terminal state, do not use a future timepoint
        %to update the Q matrix
        
        if a==2
            Q(curpos.row, curpos.col, a) = curQ + alpha*(r - curQ);
        else
            Q(curpos.row, curpos.col, a) = curQ + alpha*(r + gamma*nextQ - curQ);
        end
        
        %update time step
        curpos = nextpos; a = a_next; step = step + 1;
        
        cumReward(ei) = cumReward(ei) + r;
        
        explore_count(step) = explore;
        
%        if a_next==2
%            figure(500)
%            plot(1:ntimesteps,explore_count)
%        end



        
    end % states in each episode
    
    %Not needed
%     if (e2t_i <= length(episodesToTrack) && ei == episodesToTrack(e2t_i))
%         %add information about the environment to struct
%         episodesStruct(e2t_i).start = start;
%         episodesStruct(e2t_i).agent = agent;
%         episodesStruct(e2t_i).episode = ei;
%         episodesStruct(e2t_i).ntimesteps = ntimesteps;
%         
%         episodesStruct(e2t_i).speedbumpcol = speedbumpcol;
%         episodesStruct(e2t_i).smallreward = smallreward;
%         e2t_i = e2t_i + 1;
%     end %go to next episode in struct array
    
    
    % if the current state of the world is going to be drawn ...
%     if(selectedEpIndex <= length(plotEpisodes) && ei == plotEpisodes(selectedEpIndex))
%         curpos = start;
%         rows = []; cols = []; acts = [];
%         for i = 1:(ntimesteps) * 10,
%             %[qmax, a] = max(Q(curpos,curpos.col,:));
%             [qmax, a] = chooseAction(epsilon, Q, curpos, ntimesteps, nactions);
%             nextpos = curpos+1;
%             rows = [rows curpos];
%             cols = [cols curpos.col];
%             acts = [acts a];
% 
%             if(PosCmp(nextpos, goal) == 0), break; end
%             curpos = nextpos;
%         end % states in each episode
%         
%         %figure;
%         figure('Name',sprintf('Episode: %d', ei), 'NumberTitle','off');
%         DrawClockEpisodeState(rows, cols, acts, start, start.col, goal, goal.col, ntimesteps, fontsize, speedbumpcol, smallreward);
%         if(showTitle == 1),
%             title(sprintf('Cliff Walking %s - episode %d - (\\epsilon: %3.3f), (\\alpha = %3.4f), (\\gamma = %1.1f)', agent, ei, epsilon, alpha, gamma));
%         end
%         
%         selectedEpIndex = selectedEpIndex + 1;
%     end
    
    %plot value curves and quit counts
    if diagnos == 1
        
        figure(10002)
        subplot(4,1,1)
        plot(1:length(quits),quits) %I think this is it
        axis([0 episodeCount 0 ntimesteps]) 
        title('rt by trial')
%         hist(quits,20)
%         title('hist of quits')
        
        
        subplot(4,1,2)
        plot(smooth(Q(1:ntimesteps, 1,1)))
        hold on
        plot(smooth(Q(1:ntimesteps,1,2)),'r')
        hold off
        title('value of waits & quits')
%         
%         subplot(4,1,3)
%         plot(smooth(Q(1:ntimesteps,1, 2)))
%         title('value of quits')
%         
        subplot(4,1,3)
        scatter(nextpos.row,r)
        title('rewards gained at rt')
        hold on
        
        subplot(4,1,4)
        hist(quits,20)
        axis([0 ntimesteps 0 200])
        title('rewards gained at rt')
        
        
        %drawnow update;
    end

end % episodes loop

maxReward = -sum(cumReward);
 


function [qmax, a, explore] = chooseAction(epsilon, Q, pos, ntimesteps, nactions)
    %choose an action based on epsilon greedy strategy
    %does not permit selection of actions that move off the grid   
    wait = 1;
    quit = 2;
    
    actionNames={'wait' 'quit'};
    
    validActions = 1:nactions;

    [qmax, apos] = max(Q(pos.row, pos.col, validActions)); %qmax is the maximum value of Q, a is its index/position
    a = validActions(apos); %make sure that match above is translated into numbering for quit and wait
    explore = 0; %tracks whether current action is exploit versus explore.
%     
%     if(rand <= epsilon)
%         %explore: randomly choose one of the possible actions
%         %fprintf('Exploring\n');
%         a = validActions(randi(length(validActions), 1)); %choose randomly among valid actions
%         explore = 1;
%     end
    
% %     %Probabilistic exploration
    
    randomNormal = rand(1000,1);
     if(randomNormal(randi(length(randomNormal),1)) <= epsilon)
         explore = 1;
        if rand>(1-(1/(ntimesteps-pos.row))) %quit only 10% of the time
            a=2;
        else
            a=1;
        end
     end
    
     
     
%      if(rand<= epsilon)
%          explore = 1;
%         if rand>.9 %quit only 10% of the time
%             a=2;
%         else
%             a=1;
%         end
%     end
    
    
%    fprintf('Action chosen %s, at time point %d, explore is: %d\n', char(actionNames(a)),pos.row, explore);
