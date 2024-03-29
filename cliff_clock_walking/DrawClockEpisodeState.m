function DrawClockEpisodeState(rows, cols, acts, SRow, SCol, GRow, GCol, gridrows, gridcols, fontsize, speedbumpcol, smallreward)
% DrawCliffEpisodeState: draws a snapshot of the cliff walking world,
% specifying the sequence of actions on the world.

% written by: Sina Iravanian - June 2009
% sina@sinairv.com
% Please send your comments or bug reports to the above email address.

DrawGrid(gridrows, gridcols);
DrawTextOnCell('S', 0, SRow, SCol, gridrows, gridcols, fontsize);
DrawTextOnCell('G', 0, GRow, GCol, gridrows, gridcols, fontsize);

center = 1; %draw actions in middle of cell
color = 'k';
for i=1:length(rows),
    DrawActionOnCell(acts(i), rows(i), cols(i), gridrows, gridcols, fontsize, color, center);
end

%draw cliff -------Comment this out in new program
% for i=2:gridcols-1,
%     DrawTextOnCell('C', 0, gridrows, i, gridrows, gridcols, fontsize);
% end

%if speedbump on, then draw halfway down a given column
if speedbumpcol > 0
    for i=1:floor(0.5*gridrows)
        DrawTextOnCell('B', 0, i, speedbumpcol, gridrows, gridcols, fontsize);
    end
end

if smallreward.on == 1
    DrawTextOnCell('R', 0, smallreward.row, smallreward.col, gridrows, gridcols, fontsize);
end