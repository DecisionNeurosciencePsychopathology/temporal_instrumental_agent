%This script will create a master lookup table for each condition
%(IEV,DEV,QUADUP) the rows will be time steps the columns will be trials.
%NOTE: the TD agents sample at 100ms intervals I believe so when
%sampling occurs be sure to reduce (or repmap) the master table so it can
%sample properly. IS this inheriently a hinderence on the Q models???

ntimesteps = 500; %in 10ms bins
conds = {'IEV' 'DEV' 'QUADUP'};
ntrials = 500;

%set up seeds
global rew_rng_state;
rew_rng_seed=71;
rng(rew_rng_seed);
rew_rng_state=rng;

for i = 1:length(conds)
    %Initalize
    m.lookup = zeros(ntimesteps,ntrials);
    m.sample = ones(1,ntrials);
    for k = 1:ntrials
%         [m.lookup(:,k),~] = RewFunction(1:ntimesteps,conds{i},1);
        for j = 1:ntimesteps
            m.lookup(j,k)=RewFunction(j*10,conds{i},1); % we could do a vectorized approach but that seems to uniform?
        end
    end
    
    if i==1
        mIEV = m;
    elseif i==2
        mDEV = m;
    else
        mQUADUP = m;
    end
end

save mIEV mIEV
save mDEV mDEV
save mQUADUP mQUADUP