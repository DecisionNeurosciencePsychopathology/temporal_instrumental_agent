function b=pseudorand_rew_generator(trials, downsample)

if nargin < 2
    downsample = 10; %downsample time scale by a factor of 10 (make into centiseconds)
end

%create a trials x samples matrix so that draws are pseudorandom for each trial
nsamples = ceil(5000./downsample);
rts=0:nsamples;
rts_real=(0:nsamples).*downsample;

b.iev_rew = zeros(trials,length(rts_real));
b.dev_rew = zeros(trials,length(rts_real));
b.quadup_rew = zeros(trials,length(rts_real));

for cond = {'DEV' 'IEV', 'QUADUP'}
    for i = 1:trials
        for j = 1:length(rts_real)
            if strcmpi(cond, 'IEV')
                b.iev_rew(i,j) = RewFunction(rts_real(j),cond);
            elseif strcmpi(cond, 'DEV')
                b.dev_rew(i,j) = RewFunction(rts_real(j),cond);
            elseif strcmpi(cond, 'QUADUP')
                b.quadup_rew(i,j) = RewFunction(rts_real(j),cond);
            end
        end
    end
end

%pseudorandom uniform probabilities for sigmoid
b.randdraws=rand(1,trials);
b.rts = rts;
b.rts_real = rts_real;

end