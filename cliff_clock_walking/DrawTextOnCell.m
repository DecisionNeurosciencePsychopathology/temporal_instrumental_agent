function DrawTextOnCell(theText, rotation, row, col, gridrows, gridcols, fontsize, color, corner, jitter)
% DrawTextOnCell: draws the specified text in a specific cell of
% a grid with specified dimensions, using the rotation and fontsize
% parameters given

% written by: Sina Iravanian - June 2009
% sina@sinairv.com
% Please send your comments or bug reports to the above email address.

if nargin < 8
    color = 'k'; %default black
end

if nargin < 9
    corner = 'center';
end

if nargin < 10
    jitter = 0;
end

[xc, yc] = FindCellEdge(row, col, gridrows, gridcols, corner, jitter);


text(xc, yc, theText,  'FontSize', fontsize, 'Rotation', rotation, 'Color', color);
