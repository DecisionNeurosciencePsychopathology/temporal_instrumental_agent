function clock_options_generator(cond)
%This file creates and saves the clock_options.mat file to be used in
%clockWalker.m

%set default case to DEV
if nargin < 1
    options.cond = 'DEV';
else
    options.cond=cond;
end


options.gamma = 0.9; %discount factor: how important are are short-term versus long-term gains (0 = short-sighted, 1 = long-term payoff)
options.alpha = 0.2; %learning rate: how much to weight new information (0 = no learning, 1 = only value latest information) originally .02
options.epsilon = 0.08; %how often does agent explore? 8-10% here
options.gridcols = 50; %10ms bins so total trial length 5000ms, might have to go to full 5000 though for accurate results
options.gridrows = 2; %I need to make two options files for the degen agent and the original agent!
options.fontsize = 5;
options.showTitle = 1;

options.episodeCount = 200;
options.plotEpisodes = [20 200 700 1000 2000];
options.allowDiagonal = 0;
options.canHold = 0;

options.start.row = 1;
options.start.col = 1;
options.goal.row = 2; %I need to make two options files for the degen agent and the original agent!
options.goal.col = 10; %This probably needs to change I just wanted to get this thing running though
options.agent = 'qlearning';

options.speedbumpcol = 0; %by default, leave out speedbump
options.speedbumpcost = -10;
options.smallreward.on = 0; 
options.smallreward.row = 1;
options.smallreward.col = 1;
options.smallreward.payoff = 3;
options.snapshotEpisodes = [ 5 100 200 300 1000 2000 ]; %which episodes to track (snapshots)
options.newgoal.on = 0; %off by default
options.newgoal.episode = options.episodeCount;
options.newgoal.row = 1;
options.newgoal.col = options.gridcols; %position in NE corner
options.diag = 1; %plot diagnostics
options.decayflag =0; %incremently reduce epsilon
options.waitflag = 0; %if waiting recieves a penalty or 0

save clock_options
end
