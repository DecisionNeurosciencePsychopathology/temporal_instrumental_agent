function [x,y] = FindCellEdge(row, col, gridrows, gridcols, edge, jitter)
% FindCellCenter: finds the coordination of the center of the cell
% specified in a grid with the dimensions specified

% written by: Sina Iravanian - June 2009
% sina@sinairv.com
% Please send your comments or bug reports to the above email address.

if nargin < 6
    jitter = 0;
end

xsp = 1 / (gridcols + 2);
ysp = 1 / (gridrows + 2);
x = ((2*col + 1) / 2) * xsp;
y = 1 - (((2*row + 1) / 2) * ysp);
x = x - xsp/5;

%add/subtract half of the x and/or y spacing to place arrow between cells
switch edge
    case 'e' % east
        x = x + 0.5*xsp;
    case 's' % south
        y = y - 0.5*ysp;
    case 'w' % west
        x = x - 0.5*xsp;
    case 'n' % north
        y = y + 0.5*ysp;
    case 'ne' % northeast
        x = x + 0.5*xsp;
        y = y + 0.5*ysp;        
    case 'se' % southeast
        x = x + 0.5*xsp;
        y = y - 0.5*ysp;
    case 'sw' % southwest
        x = x - 0.5*xsp;
        y = y - 0.5*ysp;
    case 'nw' % northwest
        x = x - 0.5*xsp;
        y = y + 0.5*ysp;
    case 'h' % hold
        %fprintf('use center for a hold')        
    case 'center'
        %fprintf('plotting center')
    otherwise
        disp(sprintf('invalid action index: %d', edge))
end

%jitter using random uniform numbers within 15% of xsp or ysp
if jitter
    if any(strcmp(edge, {'n', 's', 'ne', 'nw', 'se', 'sw'}))
        x = x + .15 .*xsp .* (-1 + 2.*rand); %the code in parentheses generates a random number [-1, 1]
    end
    if any(strcmp(edge, {'e', 'w', 'ne', 'nw', 'se', 'sw'}))
        y = y + .15 .* ysp .* (-1 + 2.*rand);
    end
end