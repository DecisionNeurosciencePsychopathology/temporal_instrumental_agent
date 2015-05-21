%Make this function makeMovieModel



%Handle naming and file extension
s1 = name;
s2 = '.avi';
nameStr = [s1 s2];

%Should this be here or in the shell script, maybe call this getMovie
%instead??? I think I like that idea better do that..
writerObj = VideoWriter(nameStr); % Name it.
writerObj.FrameRate = 20; % How many frames per second.open(writerObj);
open(writerObj); %Need to open writeerObj first
for i=1:ntrials
    writeVideo(writerObj, mov(i));
end