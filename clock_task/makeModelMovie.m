%Make this function makeMovieModel
function makeModelMovie(s,agent,cond,name,reversal)


%Handle naming and file extension
s1 = name;
s2 = '.avi';
nameStr = [s1 s2];

%Only use the optimal parameters for the all condition right now
mov = getModelMovie(agent,s.(agent).opt_params(4,:),cond,reversal);

writerObj = VideoWriter(nameStr); % Name it.
writerObj.FrameRate = 20; % How many frames per second.open(writerObj);
open(writerObj); %Need to open writeerObj first
for i=1:length(mov)
    writeVideo(writerObj, mov(i));
end