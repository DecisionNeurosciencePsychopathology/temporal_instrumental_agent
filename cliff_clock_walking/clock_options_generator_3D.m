function clock_options_generator_3D(cond)
%This file creates and saves the clock_options.mat file to be used in
%clockWalker.m

%set default case to IEV
if nargin < 1
    clock_options.cond = 'IEV';
else
    clock_options.cond=cond;
end


clock_options.gamma = 0.9; %discount factor: how important are are short-term versus long-term gains (0 = short-sighted, 1 = long-term payoff)
clock_options.alpha = 0.2; %learning rate: how much to weight new information (0 = no learning, 1 = only value latest information) originally .02
clock_options.epsilon = 0.08; %how often does agent explore? 8-10% here
clock_options.gridcols = 50; %10ms bins so total trial length 5000ms, might have to go to full 5000 though for accurate results
clock_options.gridrows = 2; %I need to make two options files for the degen agent and the original agent!
clock_options.fontsize = 5;
clock_options.showTitle = 1;

clock_options.episodeCount = 200;
clock_options.plotEpisodes = [20 200 700 1000 2000];
clock_options.allowDiagonal = 0;
clock_options.canHold = 0;

clock_options.start.row = 1;
clock_options.start.col = 1;
% options.start=1;
clock_options.goal.row = 2; %I need to make two options files for the degen agent and the original agent!
clock_options.goal.col = 10; %This probably needs to change I just wanted to get this thing running though
clock_options.agent = 'qlearning';

clock_options.speedbumpcol = 0; %by default, leave out speedbump
clock_options.speedbumpcost = -10;
clock_options.smallreward.on = 0; 
clock_options.smallreward.row = 1;
clock_options.smallreward.col = 1;
clock_options.smallreward.payoff = 3;
clock_options.snapshotEpisodes = [ 5 100 200 300 1000 2000 ]; %which episodes to track (snapshots)
clock_options.newgoal.on = 0; %off by default
clock_options.newgoal.episode = clock_options.episodeCount;
clock_options.newgoal.row = 1;
clock_options.newgoal.col = clock_options.gridcols; %position in NE corner
clock_options.diag = 1; %plot diagnostics
clock_options.decayflag =0; %incremently reduce epsilon
clock_options.waitflag = 0; %if waiting recieves a penalty or 0
clock_options.ntimesteps = 50;

save clock_options
end
