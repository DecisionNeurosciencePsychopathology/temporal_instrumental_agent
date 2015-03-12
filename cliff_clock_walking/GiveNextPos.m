function nextpos = GiveNextPos(curPos, actionIndex, gridcols, gridrows)
%curPos is the agents current position (row, col)
%actionIndex is the action to take
%function moves agent from curPos to nextpos according to actionIndex
nextpos = curPos;
switch actionIndex
   case 1 % east
       nextpos.col = curPos.col + 1;
   case 2 % south
       nextpos.row = curPos.row + 1;
   case 3 % west
       nextpos.col = curPos.col - 1;
   case 4 % north
       nextpos.row = curPos.row - 1;
   case 5 % northeast 
       nextpos.col = curPos.col + 1;
       nextpos.row = curPos.row - 1;
   case 6 % southeast 
       nextpos.col = curPos.col + 1;
       nextpos.row = curPos.row + 1;
   case 7 % southwest
       nextpos.col = curPos.col - 1;
       nextpos.row = curPos.row + 1;
   case 8 % northwest
       nextpos.col = curPos.col - 1;
       nextpos.row = curPos.row - 1;
   case 9 % hold
       nextpos = curPos;
   otherwise
      disp(sprintf('invalid action index: %d', actionIndex))
end

%reset invalid moves on border of grid (results in hold action)
if(nextpos.col <= 0), nextpos.col = 1; end
if(nextpos.col > gridcols), nextpos.col = gridcols; end
if(nextpos.row <= 0), nextpos.row = 1; end
if(nextpos.row > gridrows), nextpos.row = gridrows; end
