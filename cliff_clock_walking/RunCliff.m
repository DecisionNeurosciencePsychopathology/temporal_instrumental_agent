%run cliff learning with some different setups
options.gamma = 0.9;
options.alpha = 0.05;
options.epsilon = 0.1;
options.gridcols = 10;
options.gridrows = 7;
options.fontsize = 16;
options.showTitle = 1;
options.episodeCount = 2000;
options.plotEpisodes = [20 200 700 1000 2000];
options.snapshotEpisodes = [ 5 25 100 200 300 1000 ]; %which episodes to track
options.canHold = 0;
options.start.row = 7;
options.start.col = 1;
options.goal.row = 7;
options.goal.col = 10;
options.allowDiagonal = 0;
options.speedbumpcol = 0;
options.speedbumpcost = -1;

%Q-learning -- no diagonal
options.agent = 'qlearning';
options.allowDiagonal = 0;
rq = CliffWalking(options);

figure;
%plot(1:length(rq), rq); %raw plot
plot( smooth(rq, 15, 'moving') ); %15-episode moving average

global episodesStruct;
%save episodes to movie
for i = 1:length(episodesStruct)
    epNum = episodesStruct(i).episode;
    PlotAgentMovie(episodesStruct(i), sprintf('qlearn_ep%d.avi', epNum));
end

%SARSA -- no diagonal
options.agent = 'sarsa';
options.allowDiagonal = 0;
rsarsa = CliffWalking(options);

hold on;
%figure;
%plot(1:length(rsarsa), rsarsa); %raw plot
plot( smooth(rsarsa, 15, 'moving'), 'Color', 'red' ); %15-episode moving average

global episodesStruct;
%save episodes to movie
for i = 1:length(episodesStruct)
    epNum = episodesStruct(i).episode;
    PlotAgentMovie(episodesStruct(i), sprintf('sarsa_%d.avi', epNum));
end

%Q-learning -- allow diagonal moves
options.agent = 'qlearning';
options.allowDiagonal = 1;
CliffWalking(options);

%SARSA -- allow diagonal moves
options.agent = 'sarsa';
options.allowDiagonal = 1;
CliffWalking(options);



%Q-learning -- no diagonal
%Add speedbump
options.agent = 'qlearning';
options.allowDiagonal = 0;
options.speedbumpcol = 5;
options.speedbumpcost = -5;
rq = CliffWalking(options);

global episodesStruct;
%save episodes to movie
for i = 1:length(episodesStruct)
    epNum = episodesStruct(i).episode;
    PlotAgentMovie(episodesStruct(i), sprintf('qlearn_speedbump_ep%d.avi', epNum));
end


%%%%
%small (absorbing) reward in NW corner
options.speedbumpcol = 0; %no speedbump
options.allowDiagonal = 0;
options.smallreward.on = 1;
options.smallreward.row = 1;
options.smallreward.col = 1; %northwest corner (1,1)
options.smallreward.payoff = -5; %reward

%Q-learning NW reward
options.agent = 'qlearning';
rq = CliffWalking(options);

global episodesStruct;
%save episodes to movie
for i = 1:length(episodesStruct)
    epNum = episodesStruct(i).episode;
    PlotAgentMovie(episodesStruct(i), sprintf('qlearn_nw_m5_payoff_ep%d.avi', epNum));
end

%SARSA NW reward
options.agent = 'sarsa';
rsarsa = CliffWalking(options);

global episodesStruct;
%save episodes to movie
for i = 1:length(episodesStruct)
    epNum = episodesStruct(i).episode;
    PlotAgentMovie(episodesStruct(i), sprintf('sarsa_nw_m5_payoff_ep%d.avi', epNum));
end


%%%%
%move goal from SE corner to NE corner at episode 500
options.speedbumpcol = 0;
options.smallreward.on = 0;
options.newgoal.on = 1;
options.newgoal.episode = 500;
options.newgoal.row = 1;
options.newgoal.col = options.gridcols; %position in NE corner
options.agent = 'qlearning';
options.snapshotEpisodes = [ 100 300 500 510 700 1000 ]; %which episodes to track
rqreverse = CliffWalking(options);

global episodesStruct;
%save episodes to movie
for i = 1:length(episodesStruct)
    epNum = episodesStruct(i).episode;
    PlotAgentMovie(episodesStruct(i), sprintf('qlearning_nereverse_at500_ep%d.avi', epNum));
end
