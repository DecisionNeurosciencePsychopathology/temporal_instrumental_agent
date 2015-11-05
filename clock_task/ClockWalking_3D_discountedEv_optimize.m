function [maxReward, constr, quits,cumReward,mov, ret] = ClockWalking_3D_discountedEv_optimize(options,m,rngseeds,params,plot_index,gra_options)


%Old function header 5/22/15
%[maxReward, constr, quits,cumReward,Q_vi,SARSA_vi,mov] = ClockWalking_3D_discountedEv_optimize(options,m,rngseeds,params,plot_index,gra_options)
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


%quick hack!
quits=[];
cumReward=[];



%This is for master plots
% load('s_DEV.mat')
% load('s_IEV.mat')
% load('s_QUADUP.mat')
% load('s_DEV_sarsa.mat')
% load('s_IEV_sarsa.mat')
% load('s_QUADUP_sarsa.mat')

%create folder load pwd then these, also you're going to have to load
%clock_options
current_dir = pwd;
param_path = [current_dir '\walker_params\'];


%This what cmd seems like it would be SUPER helpful in the future!
% s=what(param_path);
% matfiles=s.mat;
% for a=1:numel(matfiles)
%     load(char(matfiles(a)))
% end

% load('s_DEV_qDiscounted.mat')
% load('s_IEV_qDiscounted.mat')
% load('s_QUADUP_qDiscounted.mat')
% load('s_DEV_sarsaDiscounted.mat')
% load('s_IEV_sarsaDiscounted.mat')
% load('s_QUADUP_sarsaDiscounted.mat')

%Do you want to compare Q models to Temporal Agent?
compare_flag =1;


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
    lambda=options(4);
    
    ntimesteps = clock_options.ntimesteps;
    agent = clock_options.agent;
    episodeCount = 500; %Split 500 IEV 500 DEV or run 1000 straight
    start = clock_options.start;
    episodesToTrack = clock_options.snapshotEpisodes;
    cond = 'IEV';
    diagnos = 0;
    decay_flag = 0;
    agent='sarsa'; %Got params for Q now trying sarsa

elseif nargin ==4
    gamma = params(1);
    alpha = params(2);
    epsilon = params(3);
    lambda = params(4);
    ntimesteps = options.ntimesteps;
    episodeCount = options.episodeCount;
    start = options.start;
    agent = options.agent;
    episodesToTrack = options.snapshotEpisodes;
    smallreward = options.smallreward;
    cond = options.cond;
    diagnos = options.diag; %Reverted to set it in options
    decay_flag = options.decayflag;

    
else
    gamma = options.gamma;
    alpha = options.alpha;
    epsilon = options.epsilon;
    lambda = options.lambda; %ADD this to options
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
    diagnos = options.diag; %Whether to make plots or not
    decay_flag = options.decayflag;
    wait_punishment = options.waitflag;
    
end

%initialize movie storage
mov=repmat(struct('cdata', [], 'colormap', []), episodeCount,1);

load('mIEV.mat')
load('mDEV.mat')

% agents = {'qlearning', 'sarsa'};
% epiCount = [50, 400];
% for ii=1:length(epiCount)
%     episodeCount=epiCount(ii);
%     if ii==1
%         foo=0;
%     else
%         foo=8;
%     end
%     for grw_compare=1:length(agents)
%         agent=agents{grw_compare};
        
        
        
        
%         %Fit to optimal paramters
%         switch agent
%             case 'qlearning'
%                 if strcmp(cond,'DEV')
%                     gamma = s_DEV_qDiscounted.gamma;
%                     alpha = s_DEV_qDiscounted.alpha;
%                     epsilon = s_DEV_qDiscounted.epsilon;
%                 elseif strcmp(cond,'IEV') || strcmp(cond,'CEV') || strcmp(cond,'CEVR')
%                     gamma = s_IEV_qDiscounted.gamma;
%                     alpha = s_IEV_qDiscounted.alpha;
%                     epsilon =s_IEV_qDiscounted.epsilon;
%                 else
%                     gamma = s_QUADUP_qDiscounted.gamma;
%                     alpha = s_QUADUP_qDiscounted.alpha;
%                     epsilon = s_QUADUP_qDiscounted.epsilon;
%                 end
%             case 'sarsa'
%                 if strcmp(cond,'DEV')
%                     gamma = s_DEV_sarsaDiscounted.gamma;
%                     alpha = s_DEV_sarsaDiscounted.alpha;
%                     epsilon = s_DEV_sarsaDiscounted.epsilon;
%                 elseif strcmp(cond,'IEV') || strcmp(cond,'CEV') || strcmp(cond,'CEVR')
%                     gamma = s_IEV_sarsaDiscounted.gamma;
%                     alpha = s_IEV_sarsaDiscounted.alpha;
%                     epsilon = s_IEV_sarsaDiscounted.epsilon;
%                 else
%                     gamma = s_QUADUP_sarsaDiscounted.gamma;
%                     alpha = s_QUADUP_sarsaDiscounted.alpha;
%                     epsilon = s_QUADUP_sarsaDiscounted.epsilon;
%                 end
%         end
%         


        %Load in look up tables
%         if strcmp(cond,'DEV')
%            load('mDEV.mat');
%            m = mDEV;
%         elseif strcmp(cond,'IEV') || strcmp(cond,'CEV') || strcmp(cond,'CEVR')
%            load('mIEV.mat');
%            m = mIEV;
%         else
%            load('mQUADUP.mat');
%            m = mQUADUP;
%         end


        
        %states for random generators are shared across functions to allow for repeatable draws
        global rew_rng_state explore_rng_state;
        %rngseeds=[98 83]; %changed this up
        
        
        %initialize states for two repeatable random number generators using different seeds
        rew_rng_seed=rngseeds(1);
        explore_rng_seed=rngseeds(2);
        rng(rew_rng_seed);
        rew_rng_state=rng;
        rng(explore_rng_seed);
        explore_rng_state=rng;
        
        fprintf('running agent with gamma: %.3f, alpha: %.3f, epsilon: %.3f, lambda: %.3f and rngseeds: %s \n', gamma, alpha, epsilon, lambda, num2str(rngseeds))
        
        %The amount of decay you want to decrease epsilon by
        decay_val=.95; %Currently 5%
        
        global four_panel;
        four_panel = 1;
        
        
        selectedEpIndex = 1;
        %quit or wait
        nactions = 2;
        
        % initialize Q with zeros (no initial information about environment)
        % Q is the action-value function: expected reward for each possible action at a given timestep in the trial.
        global Q;
        Q = zeros(ntimesteps, nactions, nactions);
        
        %Set up eligibility traces
        e = zeros(ntimesteps,ntimesteps);
        for t = 1:ntimesteps
            for w = t:ntimesteps
                e(w,t)=(gamma^(w-t))*(lambda^(w-t));
            end
        end
        
        
        
        %What if Q was intially 1
        %Q = ones(ntimesteps, nactions, nactions);
        
        ah = ones(ntimesteps,1);
        cumReward = zeros(episodeCount, 1); %vector to track total reward earned per episode
        
        a = 0; % an invalid action
        % loop through episodes
        % agent is allowed to sample environment many times and carry forward its learning to subsequent episodes.
        
        %To show the updates in Q as the agent moves through environment, take full snapshots of some episodes
        e2t_i = 1; %index of episode being tracked
        episodesStruct = []; %will be a vector (episodes) of structs containing cell arrays of Q and a (per action/step)
        
        %set factor depending upon a 10ms or 100ms bin size
        factor=100; %bins are 100ms
        
        quits=nan(episodeCount,1);
        
        explore_count=zeros(ntimesteps,1);
        
        %Initalize dot sizes based on rts
        dot_size = 30*(ones(length(quits),2));
        
        %Reversal hack for switching from IEV to DEV
        rev_go=options.reversal_go;
        
        %Expected value
        ev_i = zeros(1, episodeCount);
        
        
        for ei = 1:episodeCount,
            
            %Reversal hack
%             if (rev_go ==1) && ei == episodeCount/2+1
%                 vperm_run = m.vperm_run; %Currently this only exsists for the reversal runs
%                 if strcmp(m.name, 'IEV')
%                     m=mDEV; %if it is IEV after x trials switch
%                     m.lookup = m.lookup(:,vperm_run);
%                     cond = m.name;
%                 else
%                     m=mIEV; %else it is DEV after x trials switch to IEV
%                     m.lookup = m.lookup(:,vperm_run);
%                     cond = m.name;
%                 end
%             end
            
            
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
                ah(step) = a;
                
                
                % 2) take action a to move from curpos to nextpos
                % 3) receive reward r (at nextpos)
                nextpos.row = curpos.row+1; %next timestep
                %fprintf('ei is: %d, nextpos is %d\n', ei, nextpos);
                %fprintf('Location: %d %d\n', nextpos, nextpos);
                if a == 2 %i.e., quit
                    %if we have arrived at the goal state, mark as finished and provide reward 0
                    episodeFinished = 1;
                    %r = RewFunction(nextpos.row.*factor, cond);
                    [r, ev_i(ei), m] = getNextRew(nextpos.row.*10, m);
%                     [r, m] = getNextRew(nextpos.row.*10, m);
%                     [~,ev_i(ei)] = RewFunction(nextpos.row.*factor, cond);
                    
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
                    ah(step+1)=a_next;
                end
                
                % update Q matrix (action-value)
                % learned value is [r + gamma*nextQ]
                % old value is curQ
                
                %I think this is what Alex wanted, once the agent has decided it
                %wants to go into the terminal state, do not use a future timepoint
                %to update the Q matrix
                
                %set t to curpos.row
                t = curpos.row;
                
                
                %For SARSA
                
                
                if a==2
                    delta = r-curQ;
                    Q(t,curpos.col,a) = Q(t,curpos.col,a) + alpha*delta;
                    if t==1
                        %Q(:,curpos.col,a-1) = Q(:,curpos.col,a-1) + alpha*delta*e(t,:)'; %This was the line i was not sure about being correct
                    else
                        Q(:,curpos.col,a-1) = Q(:,curpos.col,a-1) + alpha*delta*e(t-1,:)';
                    end
                else
                    delta = r+gamma*nextQ-curQ;
                    Q(:,curpos.col,a) = Q(:,curpos.col,a) + alpha*delta*e(t,:)';
                end
                
                %See where Q quit is greater than Q wait
                Q_quit_grt_wait(ei,:) = Q(:,1,1)<Q(:,1,2);
                
                %
                %         if a==2
                %             nextQ=0;
                %         end
                %             delta = r+gamma*nextQ-curQ;
                %             Q(:,curpos.col,a) = Q(:,curpos.col,a) + alpha*delta*e(t,:)';
                
                
                
                %         for w = 1:t
                %             e(ei,w) = (gamma^(t-w))*(lambda^(t-w));
                %             curQ = Q(w, curpos.col, ah(w));
                %             nextQ = Q(w+1, curpos.col, ah(w+1));
                %             if a==2 && w==t
                %                 Q(w, curpos.col, ah(w)) = curQ + alpha*(r - curQ);
                %                 %break %Do we terminate right after a quit?
                %             else
                %                 Q(w, curpos.col, ah(w)) = curQ + alpha*(r + gamma*nextQ - curQ)*e(ei,w);
                %             end
                %         end
                
                
                
                %update time step
                curpos = nextpos; a = a_next; step = step + 1;
                
                cumReward(ei) = cumReward(ei) + r;
                
                quit_rews(ei) = r;
                
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
                
                figure(71)
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
                subplot(4,1,4)
                scatter(curpos.row,r)
                title('rewards gained at rt')
                hold on
                
                %                 subplot(4,1,4)
                %                 hist(quits,20)
                %                 axis([0 ntimesteps 0 200])
                %                 title('rewards gained at rt')
                
                subplot(4,1,3)
                %plot(Q(1:ntimesteps,1,2)-Q(1:ntimesteps, 1,1))
                %plot(smooth(Q_quit_grt_wait))
                plot(e(t,:))
                %hold on
                %axis([0 ntimesteps 0 200])
                title('eligibility trace')
                
                mov(ei) = getframe(gcf);
                
            elseif diagnos == 2
                %plots q-learnfing first then sarsa
                
                quit_idx = zeros(length(quits),1);
                
                %%%This info is used to create subplots
                rows=gra_options.rows;
                cols=gra_options.cols;
                font_size = gra_options.font_size;
                line_width = gra_options.line_width;
                fig_num = gra_options.fig_numbers;
                bfor_rev = gra_options.bfor_rev;
                add_row=0;
                
                if bfor_rev==1 && add_row==0
                    rows=rows+1; %add an additional middle row show value/uncertainity/choice bfore reversal
                    add_row=1;
                end
                %%%%%%%%%
                
                %Split graph into seperate panels
                if foo<3
                    figure(fig_num(1))
                else
                    figure(fig_num(2))
                end
                set(gca,'FontSize',18);
                
                if strcmp(agent,'qlearning')
                    
                    %size=14;
                    subplot(rows,cols,plot_index(1))
                    plot(1:length(quits),quits*10,'Linewidth',2) %I think this is it
                    axis([0 episodeCount 0 ntimesteps*10+50])
                    %ylabel('Time x 10ms','FontSize',size)
                    %xlabel('Trials','FontSize',size)
                    %title('RT by Trial','FontWeight','bold','FontSize',size)
                    
                    
                    %See this is why you need to follow the DRY
                    %principle!!!!
                    if bfor_rev==1
                        subplot(rows,cols,plot_index(1)+cols)
                        if ei<episodeCount/2
                            plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps, 1,1)),'r','Linewidth',line_width)
                            hold on
                            plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps,1,2)),'Linewidth',line_width)
                            hold on
                            quit_idx=~isnan(quits(1:ei));
                            scatter(quits(quit_idx)*10,quit_rews(quit_idx),dot_size(quit_idx,2),'k','filled')
                            if ~isnan(quits(ei))
                                dot_size(quits(ei),2) = dot_size(quits(ei),2)+2;
                            end
                            hold off
                        end
                    end
                    
                    
                    subplot(rows,cols,plot_index(1)+cols*(rows-1));
                    %plotyy(1:50,smooth(Q(1:ntimesteps, 1,1)),1:50,smooth(Q(1:ntimesteps,1,2)))
                    %hold on
                    plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps, 1,1)),'r','Linewidth',line_width)
                    hold on
                    plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps,1,2)),'Linewidth',line_width)
                    hold on
                    
                    quit_idx=~isnan(quits(1:ei));
                    scatter(quits(quit_idx)*10,quit_rews(quit_idx),dot_size(quit_idx,1),'k','filled')
                    if ~isnan(quits(ei))
                        dot_size(quits(ei),1) = dot_size(quits(ei),1)+2;
                    end
                    
                    hold off
                    %ylabel('Reward','FontSize',size)
                    %xlabel('Time x 10ms','FontSize',size)
                    %title(sprintf('Q_w_a_i_t(blue) Q_q_u_i_t(red)\n Rewards Sampled(circles)'),'FontWeight','bold','FontSize',size)
                    if foo>0
                        hh=sprintf('Q-Learning');
                        %hh=sprintf('Q-learning: Alpha %.3f, Gamma %.3f,\n Epsilon %.3f, Lambda %.3f',alpha,...
                        %    gamma, epsilon, lambda);
                        %xlabel({'Time x 10ms';hh},'FontSize',size)
                        xlabel({hh},'FontSize',font_size);
                    end
                    Q_vi(ei,:) = Q(1:ntimesteps,1,2);
                else
                    %size=14;
                    subplot(rows,cols,plot_index(2))
                    plot(1:length(quits),quits*10,'Linewidth',2) %I think this is it
                    axis([0 episodeCount 0 ntimesteps*10+50])
                    %ylabel('Time x 10ms','FontSize',size)
                    %xlabel('Trials','FontSize',size)
                    %title('RT by Trial','FontWeight','bold','FontSize',size)
                    
                    
                    if bfor_rev==1
                        subplot(rows,cols,plot_index(2)+cols)
                        if ei<episodeCount/2
                            plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps, 1,1)),'r','Linewidth',line_width)
                            hold on
                            plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps,1,2)),'Linewidth',line_width)
                            hold on
                            quit_idx=~isnan(quits(1:ei));
                            scatter(quits(quit_idx)*10,quit_rews(quit_idx),dot_size(quit_idx,2),'k','filled')
                            if ~isnan(quits(ei))
                                dot_size(quits(ei),2) = dot_size(quits(ei),2)+2;
                            end
                            hold off
                        end
                    end
                    
                    
                    subplot(rows,cols,plot_index(2)+cols*(rows-1));
                    %plotyy(1:50,smooth(Q(1:ntimesteps, 1,1)),1:50,smooth(Q(1:ntimesteps,1,2)))
                    %hold on
                    plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps, 1,1)),'r','Linewidth',line_width)
                    hold on
                    plot((1:ntimesteps)*10,smooth(Q(1:ntimesteps,1,2)),'Linewidth',line_width)
                    hold on
                    quit_idx=~isnan(quits(1:ei));
                    scatter(quits(quit_idx)*10,quit_rews(quit_idx),dot_size(quit_idx,1),'k','filled')
                    if ~isnan(quits(ei))
                        dot_size(quits(ei),1) = dot_size(quits(ei),1)+2;
                    end
                    hold off
                    %ylabel('Reward','FontSize',size)
                    %xlabel('Time x 10ms','FontSize',size)
                    %title(sprintf('Q_w_a_i_t(blue) Q_q_u_i_t(red)\n Rewards Sampled(circles)'),'FontWeight','bold','FontSize',size)
                    if foo>0
                        hh=sprintf('SARSA');
                        %                         hh=sprintf('SARSA: Alpha %.3f, Gamma %.3f,\n Epsilon %.3f, Lambda %.3f',alpha,...
                        %                             gamma, epsilon, lambda);
                        %xlabel({'Time x 10ms';hh},'FontSize',size)
                        xlabel({hh},'FontSize',font_size);
                    end
                    SARSA_vi(ei,:) = Q(1:ntimesteps,1,2);
                end
                
                %                 clock_logistic_operator(params,[93 83],cond,ntimesteps)
                
                
                %drawnow update;
                %Update dot size
                %dot_size(quits(ei)) = dot_size(quits(ei))+2;
            end
            
        end % episodes loop
        
        maxReward = -sum(cumReward);
        %update plot index
        
        ret.ev_i = ev_i;
        ret.rts = quits; 
        ret.rew_i = cumReward;
        
    end
    %     plot_index(1)=plot_index(1)+(cols*2);
    %     plot_index(2)=plot_index(2)+(cols*2);
    
% end
%TA_opt = [-0.045 0.2 .9706];
% TA_opt = [-0.17 0.2 .01]; %Optimal params w/wo noise
% % %TA_kal_opt = [.055 .43]; %Via simps multi run
% % TA_kal_opt = [.479 .0363];
% % TA_kal_opt_wExploreNoise = [0.42161 0.041876];
% % TA_kal_opt_BRUTAL = [.17 .17];
% % TA_kal_opt_GA = [.4764 .0305];
% 
% kalman_params = [.6 .05 .5 .2];
% kalman_seeds = [69 56 20 30];
% 
% seeds = [99 88];
% seeds_GA = [280 432];
% if strcmp(cond,'IEV') && compare_flag==1
%     %load('s_IEV_TA.mat')
%     %clock_logistic_operator(s_IEV_TA.parameters,[98 83],'IEV')
%     clock_logistic_operator_plots(TA_opt,[98 83],'IEV')
% elseif strcmp(cond,'DEV')
%     %load('s_DEV_TA.mat')
%     %clock_logistic_operator(s_DEV_TA.parameters,[98 83],'DEV')
%     clock_logistic_operator_plots(TA_opt,[98 83],'DEV')
% else
%     %load('s_QUADUP_TA.mat')
%     %clock_logistic_operator(s_QUADUP_TA.parameters,[98 83],'QUADUP')
%    clock_logistic_operator_plots(TA_opt,[98 83],'QUADUP')
% end

%Compare selective and random walk exploration versions of kalman filter


%Old code
% for grw_compare=1:2
%     if grw_compare==1
%         kalman_params=[.6 .05 .9 .2]; %highly GRW
%     else
%         kalman_params=[.6 .05 0 .2]; %solely stratego
%     end
%     if strcmp(cond,'IEV') && compare_flag==1
%         %Same seeds
%         clock_logistic_operator_kalman_plotter(kalman_params,[339 295 20 30],'IEV',80, 24, 500, 2,grw_compare)
%         %new seeds
%         %clock_logistic_operator_kalman_plotter([.370 .21381],[430 47],'IEV')
%     elseif strcmp(cond,'DEV')
%         %Same seeds
%         clock_logistic_operator_kalman_plotter(kalman_params,[471 376 20 30],'DEV',80, 24, 500, 2,grw_compare)
%         %new seeds
%         %clock_logistic_operator_kalman_plotter([.370 .21381],[430 47],'DEV')
%     elseif strcmp(cond,'CEV')
%         %Same seeds
%         clock_logistic_operator_kalman_plotter(kalman_params,kalman_seeds,'CEV',80, 24, 500, 2,grw_compare)
%     elseif strcmp(cond,'CEVR')
%         %Same seeds
%         clock_logistic_operator_kalman_plotter(kalman_params,kalman_seeds,'CEVR',80, 24, 500, 2,grw_compare)
%     else
%         %Same seeds
%         clock_logistic_operator_kalman_plotter(kalman_params,[467 156 20 30],'QUADUP',80, 24, 500, 2,grw_compare)
%         %new seeds
%         %clock_logistic_operator_kalman_plotter([.370 .21381],[430 47],'QUADUP')
%     end
% end



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
    if rand>(1-(1/(ntimesteps-pos.row)))
        a=2;
    else
        a=1;
    end
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
