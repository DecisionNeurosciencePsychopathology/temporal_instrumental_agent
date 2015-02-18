function DrawAgentOnCell(row, col, gridrows, gridcols, fontsize, color)
% DrawActionOnCell: draws the symbol of different grid-world actions on the
% cell specified

% written by: Sina Iravanian - June 2009
% sina@sinairv.com
% Please send your comments or bug reports to the above email address.

rotation = 0;
textToDraw = 'o';

if nargin < 6
    color = 'r';
end

DrawTextOnCell(textToDraw, rotation,  row, col, gridrows, gridcols, fontsize, color);
