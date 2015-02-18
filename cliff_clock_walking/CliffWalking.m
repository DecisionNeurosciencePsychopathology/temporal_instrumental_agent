function [cumReward] = CliffWalking(options)
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

if(nargin < 1),
    gamma = 0.9; %discount factor: how important are are short-term versus long-term gains (0 = short-sighted, 1 = long-term payoff)
    alpha = 0.05; %learning rate: how much to weight new information (0 = no learning, 1 = only value latest information)
    epsilon = 0.1; %how often does agent explore? 10% here
    gridcols = 10;
    gridrows = 7;
    fontsize = 16;
    showTitle = 1;

    episodeCount = 2000;
    plotEpisodes = [20 200 700 1000 2000];
    allowDiagonal = 0; 
    canHold = 0;

    start.row = 7;
    start.col = 1;
    goal.row = 7;
    goal.col = 10;
    agent = 'qlearning';
    speedbumpcol = 0; %by default, leave out speedbump
    speedbumpcost = -10;
    smallreward.on = 0;
    smallreward.row = 1;
    smallreward.col = 1;
    smallreward.payoff = 3;
    episodesToTrack = [ 5 100 200 300 1000 2000 ]; %which episodes to track (snapshots)
    newgoal.on = 0; %off by default
    newgoal.episode = episodeCount;
    newgoal.row = 1;
    newgoal.col = gridcols; %position in NE corner
else
    gamma = options.gamma;
    alpha = options.alpha;
    epsilon = options.epsilon;
    gridcols = options.gridcols; 
    gridrows = options.gridrows;
    fontsize = options.fontsize;
    showTitle = options.showTitle;

    episodeCount = options.episodeCount;
    plotEpisodes = options.plotEpisodes;
    allowDiagonal = options.allowDiagonal; 
    canHold = options.canHold;

    start = options.start;
    goal = options.goal;
    agent = options.agent;
    speedbumpcol = options.speedbumpcol;
    speedbumpcost = options.speedbumpcost;
    episodesToTrack = options.snapshotEpisodes;
    smallreward = options.smallreward;
    newgoal = options.newgoal;
end

selectedEpIndex = 1;
if(allowDiagonal ~= 0),  actionCount = 8; else actionCount = 4; end
if(canHold ~= 0 && allowDiagonal ~= 0), actionCount = actionCount + 1; end

% initialize Q with zeros (no initial information about environment)
% Q is the action-value function: expected reward for each possible action at a given position in the maze.
global Q;
Q = zeros(gridrows, gridcols, actionCount);
cumReward = zeros(episodeCount, 1); %vector to track total reward earned per episode

a = 0; % an invalid action
% loop through episodes
% agent is allowed to sample environment many times and carry forward its learning to subsequent episodes.

%To show the updates in Q as the agent moves through environment, take full snapshots of some episodes
e2t_i = 1; %index of episode being tracked
episodesStruct = []; %will be a vector (episodes) of structs containing cell arrays of Q and a (per action/step)

reversalOccurred = 0; %whether goal has been moved

for ei = 1:episodeCount,
    fprintf('Running episode %d\n', ei);
    curpos = start; %initialize current position to the start
    nextpos = start; %initialize the next position to the start
    
    % 1) choose initial action at current position using an epsilon greedy policy derived from Q
    [qmax, a, explore] = chooseAction(epsilon, Q, curpos, gridrows, gridcols, actionCount);
    
    %switch goal to newgoal if requested
    if reversalOccurred == 0 && newgoal.on == 1 && ei >= newgoal.episode 
        goal = newgoal;
        reversalOccurred = 1;
    end
    
    episodeFinished = 0;
    %Wander in the environment until agent 1) achieves the goal state, or 2) falls off the cliff.
    %This constitutes a single episode (i.e., arrive at the final/absorbing state).
        
    step = 1;
    while(episodeFinished == 0)
        % 2) take action a to move from curpos to nextpos
        % 3) receive reward r (at nextpos)
        nextpos = GiveNextPos(curpos, a, gridcols, gridrows);
        %fprintf('Location: %d %d\n', nextpos.row, nextpos.col);
        if(PosCmp(nextpos, goal) == 0)
            %if we have arrived at the goal state, mark as finished and provide reward 0
            episodeFinished = 1;
            r = 0;
        elseif(speedbumpcol > 0 && nextpos.col == speedbumpcol && nextpos.row <= floor(0.5 .* gridrows))
            %put a speedbump on the upper half of the rows at a given column (speedbumpcol)
            %this is a non-absorbing negative reward (to be avoided)
            r = speedbumpcost;
        elseif (smallreward.on == 1 && PosCmp(nextpos, smallreward) == 0)
            %provide agent a small reward (absorbing)
            episodeFinished = 1;
            r = smallreward.payoff;
        elseif(nextpos.row == gridrows && 1 < nextpos.col && nextpos.col < gridcols)
            %if the next move puts us on the bottom row and is 1 < column < gridcols,
            %then we have fallen off the cliff (bottom row, middle columns)
            episodeFinished = 1;
            r = -100; %provide a big punishment
        else
            r = -1; %every move costs 1
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
        [qmax, a_next, explore] = chooseAction(epsilon, Q, nextpos, gridrows, gridcols, actionCount);
        
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
        Q(curpos.row, curpos.col, a) = curQ + alpha*(r + gamma*nextQ - curQ);    
        
        %update time step
        curpos = nextpos; a = a_next; step = step + 1;
        
        cumReward(ei) = cumReward(ei) + r;
               
    end % states in each episode
    
    if (e2t_i <= length(episodesToTrack) && ei == episodesToTrack(e2t_i))
        %add information about the environment to struct
        episodesStruct(e2t_i).start = start;
        episodesStruct(e2t_i).goal = goal;
        episodesStruct(e2t_i).agent = agent;
        episodesStruct(e2t_i).episode = ei;
        episodesStruct(e2t_i).gridrows = gridrows;
        episodesStruct(e2t_i).gridcols = gridcols;
        episodesStruct(e2t_i).speedbumpcol = speedbumpcol;
        episodesStruct(e2t_i).smallreward = smallreward;
        e2t_i = e2t_i + 1;
    end %go to next episode in struct array
    
    % if the current state of the world is going to be drawn ...
    if(selectedEpIndex <= length(plotEpisodes) && ei == plotEpisodes(selectedEpIndex))
        curpos = start;
        rows = []; cols = []; acts = [];
        for i = 1:(gridrows + gridcols) * 10,
            %[qmax, a] = max(Q(curpos.row,curpos.col,:));
            [qmax, a] = chooseAction(epsilon, Q, curpos, gridrows, gridcols, actionCount);
            nextpos = GiveNextPos(curpos, a, gridcols, gridrows);
            rows = [rows curpos.row];
            cols = [cols curpos.col];
            acts = [acts a];

            if(PosCmp(nextpos, goal) == 0), break; end
            curpos = nextpos;
        end % states in each episode
        
        %figure;
        figure('Name',sprintf('Episode: %d', ei), 'NumberTitle','off');
        DrawCliffEpisodeState(rows, cols, acts, start.row, start.col, goal.row, goal.col, gridrows, gridcols, fontsize, speedbumpcol, smallreward);
        if(showTitle == 1),
            title(sprintf('Cliff Walking %s - episode %d - (\\epsilon: %3.3f), (\\alpha = %3.4f), (\\gamma = %1.1f)', agent, ei, epsilon, alpha, gamma));
        end
        
        selectedEpIndex = selectedEpIndex + 1;
    end

end % episodes loop

function c = PosCmp(pos1, pos2)
%are the two positions identical? (i.e., 0)
c = pos1.row - pos2.row;
if(c == 0)
    c = c + pos1.col - pos2.col;
end

% deprecated, just use randi function directly
% function n = IntRand(lowerBound, upperBound)
%     %goal is to generate random uniform integers
%     %note that as previously coded, would only generate numbers [lowerBound, upperBound - 1]
%     upperBound = upperBound + 1; %ensure that numbers range [lowerBound, upperBound]
%     n = floor((upperBound - lowerBound) * rand + lowerBound);

function [qmax, a, explore] = chooseAction(epsilon, Q, pos, gridrows, gridcols, actionCount)
    %choose an action based on epsilon greedy strategy
    %does not permit selection of actions that move off the grid   
    east = 1;
    south = 2;
    west = 3;
    north = 4;
    northeast = 5;
    southeast = 6;
    southwest = 7;
    northwest = 8;
    hold = 9;
    
    actionNames={'east' 'south' 'west' 'north' 'northeast' 'southeast' 'southwest' 'northwest' 'hold'};
    
    validActions = 1:actionCount;
    
    %restrict max Q selection to valid moves
    % gridrows are numbered top-to-bottom, gridcols are numbered left-to-right.
    if pos.row == gridrows, validActions = validActions(validActions ~= south); end %can't go south on bottom row
    if pos.row == 1, validActions = validActions(validActions ~= north); end %can't go north on top row
    if pos.col == gridcols, validActions = validActions(validActions ~= east); end %can't go east on right col
    if pos.col == 1, validActions = validActions(validActions ~= west); end %can't go west on left col

    if actionCount > 4
        if pos.row == gridrows, validActions = validActions(~ismember(validActions, [southeast southwest])); end %can't go SE or SW on bottom row
        if pos.row == 1, validActions = validActions(~ismember(validActions, [northeast northwest])); end %can't go NE or NW on top row
        if pos.col == gridcols, validActions(~ismember(validActions, [southeast northeast])); end %can't go NE or SE on right col
        if pos.col == 1, validActions(~ismember(validActions, [southwest northwest])); end %can't go NW or SW on left col
    end  
    
    [qmax, apos] = max(Q(pos.row, pos.col, validActions)); %qmax is the maximum value of Q, a is its index/position
    a = validActions(apos); %make sure that match above is translated into the n, e, s, w nomenclature above
    explore = 0; %tracks whether current action is exploit versus explore.
    
    if(rand <= epsilon)
        %explore: randomly choose one of the possible actions
        %fprintf('Exploring\n');
        a = validActions(randi(length(validActions), 1)); %choose randomly among valid actions
        explore = 1;
    end
    
    %fprintf('Action chosen %s\n', char(actionNames(a)));