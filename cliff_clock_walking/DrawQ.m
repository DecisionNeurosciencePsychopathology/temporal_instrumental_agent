function DrawQ(Q, row, col, gridrows, gridcols, fontsize, color)

%draw Q matrix overlay around agent
if nargin < 6, fontsize = 11; end
if nargin < 7, color = 'r'; end

if size(Q,3) == 4,
    %n, e, s, w
    Qeast = Q(row, col, 1);
    Qsouth = Q(row, col, 2);
    Qwest = Q(row, col, 3);
    Qnorth = Q(row, col, 4);
    
    DrawTextOnCell(num2str(Qeast, 1), 0, row, col+1, gridrows, gridcols, fontsize, color, 'ne');
    DrawTextOnCell(num2str(Qsouth, 1), 0, row+1, col, gridrows, gridcols, fontsize, color, 'ne');
    DrawTextOnCell(num2str(Qwest, 1), 0, row, col-1, gridrows, gridcols, fontsize, color, 'ne');
    DrawTextOnCell(num2str(Qnorth, 1), 0, row-1, col, gridrows, gridcols, fontsize, color, 'ne');
end

