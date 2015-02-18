function DrawActionOnCell(actionIndex, row, col, gridrows, gridcols, fontsize, color, center, jitter)
% DrawActionOnCell: draws the symbol of different grid-world actions on the
% cell specified

% written by: Sina Iravanian - June 2009
% sina@sinairv.com
% Please send your comments or bug reports to the above email address.

rotation = 0;
textToDraw = 'o';

if nargin < 7
    color = 'k';
end

if nargin < 8
    center = 1 ; %default to draw in center
end

if nargin < 9
    jitter = 0 ; %default to no jitter
end


switch actionIndex
   case 1 % east
       textToDraw = '\rightarrow';
       rotation = 0;
       corner = 'e';
   case 2 % south
       textToDraw = '\downarrow';
       rotation = 0;
       corner = 's';
   case 3 % west
       textToDraw = '\leftarrow';
       rotation = 0;
       corner = 'w';
   case 4 % north
       textToDraw = '\uparrow';
       rotation = 0;
       corner = 'n';
   case 5 % northeast 
       textToDraw = '\rightarrow';
       rotation = 45;
       corner = 'ne';
   case 6 % southeast 
       textToDraw = '\downarrow';
       rotation = 45;
       corner = 'se';
   case 7 % southwest
       textToDraw = '\leftarrow';
       rotation = 45;
       corner = 'sw';
   case 8 % northwest
       textToDraw = '\uparrow';
       rotation = 45;
       corner = 'nw';
   case 9 % hold
       textToDraw = 'o';
       rotation = 0;
       corner = 'h';
   otherwise
      disp(sprintf('invalid action index: %d', actionIndex))
end

if center == 1
    corner = 'center'; %tell drawing routine to center
end

DrawTextOnCell(textToDraw, rotation,  row, col, gridrows, gridcols, fontsize, color, corner, jitter);
