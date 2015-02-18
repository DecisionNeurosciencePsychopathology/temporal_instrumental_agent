function PlotAgentMovie(ep, fname, fontsize)
%show the agent's estimate of Q as it moves through the environment for a given episode
%ep is a struct containing two cell arrays: Q and a.
%the length of these arrays is the number of actions/steps in the episode

if nargin < 3
    fontsize = 16;
end

nsteps=length(ep.Q);
figure('Name',sprintf('Episode: %d', ep.episode), 'NumberTitle','off');

DrawGrid(ep.gridrows, ep.gridcols); %grid
DrawTextOnCell('S', 0, ep.start.row, ep.start.col, ep.gridrows, ep.gridcols, fontsize); %start
DrawTextOnCell('G', 0, ep.goal.row, ep.goal.col, ep.gridrows, ep.gridcols, fontsize); %goal

%draw cliff
for i=2:ep.gridcols-1,
    DrawTextOnCell('C', 0, ep.gridrows, i, ep.gridrows, ep.gridcols, fontsize);
end

%if speedbump on, then draw halfway down a given column
if ep.speedbumpcol > 0
    for i=1:floor(0.5*ep.gridrows)
        DrawTextOnCell('B', 0, i, ep.speedbumpcol, ep.gridrows, ep.gridcols, fontsize);
    end
end

%if small reward on, plot its location
if ep.smallreward.on == 1
    DrawTextOnCell('R', 0, ep.smallreward.row, ep.smallreward.col, ep.gridrows, ep.gridcols, fontsize);
end

%preallocate structure
mov(1:nsteps) = struct('cdata', [],...
                        'colormap', []);

%plot agent moving through environment (play actions)
%only show 
curpos=ep.start;
center = 0; %show actions on boundaries between cells
jitter = 1; %jitter arrows
for i=1:nsteps,
    %get Q, a, and explore at this step
    a = ep.a{i};
    explore = ep.explore{i};
    Q = ep.Q{i};
    
    if explore == 1
        color='b'; %draw explore actions as blue
    else
        color='k';
    end
    
    %obtain a figure snapshot with and without Q overlay
    
    
    %to obtain row and col positions, play out actions
    nextpos = GiveNextPos(curpos, a, ep.gridcols, ep.gridrows);
    
    DrawActionOnCell(a, curpos.row, curpos.col, ep.gridrows, ep.gridcols, fontsize, color, center, jitter);
    DrawAgentOnCell(curpos.row, curpos.col, ep.gridrows, ep.gridcols, fontsize, color);
    
    %beforeQ = copyobj(gcf, 0);
    
    %leave off Q for now
    %DrawQ(Q, curpos.row, curpos.col, ep.gridrows, ep.gridcols);
    
    mov(i) = getframe(gcf);
    curpos = nextpos; %increment action step
    %clf;
    %beforeQ;
end
%movie(mov, 10);
movie2avi(mov, fname, 'compression', 'None', 'FPS', 5);