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


function mov = getModelMovie(agent, params,condition,reversal)

%load in data files
load('mIEV.mat')
load('mDEV.mat')
load('mQUADUP.mat')
load('clock_options.mat')
load('s.mat')

clock_options.diag=1;

if strcmp(condition,'IEV')
    m = mIEV;
elseif strcmp(condition, 'DEV')
    m = mDEV;
else
    m = mQUADUP;
end

%Permute the reward matrix
runseed = 55;
rng(runseed);
%seeds=randi([1 500], 1, 4);
seeds = [98 83 66 10];
%vperm = zeros(1,length(mIEV.sample));
%vperm = randperm(length(mIEV.lookup));
vperm_run = randperm(length(mIEV.sample));
m.vperm_run = vperm_run;
m.lookup = m.lookup(:,vperm_run); %permute here as well
%m.lookup = m.lookup(vperm,:); %Permute

%Trial set up variables
ntrials = 200;
nbasis = 24;
ntimesteps = 500;

%Reversal
clock_options.reversal_go=reversal;
rev = reversal;

switch agent
    case 'kalmanUV'
        [~,~,~,mov,~] = clock_logistic_operator_kalman_optimize(params, seeds, condition, ntrials, nbasis, ntimesteps, 1, 1, m,rev);
    case 'kalmanLogistic'
        [~,~,~,mov,~] = clock_logistic_operator_kalman_optimize([params 0 0.2], seeds, condition, ntrials, nbasis, ntimesteps, 1, 0, m,rev);
    case 'kalmanGRW'
        [~,~,~,mov,~] = clock_logistic_operator_kalman_optimize([params 0.9 0.2], seeds, condition, ntrials, nbasis, ntimesteps, 1, 0, m,rev);
    case 'qlearning'
        clock_options.agent = 'qlearning';
       [~,~,~,~,mov] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(1:2),params);
    case 'sarsa'
        clock_options.agent = 'sarsa';
        [~,~,~,~,mov] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(1:2),params);
end