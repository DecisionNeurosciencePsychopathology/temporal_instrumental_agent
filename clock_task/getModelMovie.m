function mov = getModelMovie(agent, params,condition)

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
seeds=randi([1 500], 1, 4);
vperm = zeros(1,length(mIEV.sample));
vperm = randperm(length(mIEV.lookup));
m.lookup = m.lookup(vperm,:); %Permute

%Trial set up variables
ntrials = 200;
nbasis = 24;
ntimesteps = 500;


switch agent
    case 'kalmanUV'
        [~,~,~,mov,~] = clock_logistic_operator_kalman_optimize(params, seeds, condition, ntrials, nbasis, ntimesteps, 1, 1, m);
    case 'kalmanLogistic'
        [~,~,~,mov,~] = clock_logistic_operator_kalman_optimize([params 0 0.2], seeds, condition, ntrials, nbasis, ntimesteps, 1, 0, m);
    case 'kalmanGRW'
        [~,~,~,mov,~] = clock_logistic_operator_kalman_optimize([params 0.9 0.2], seeds, condition, ntrials, nbasis, ntimesteps, 1, 0, m);
    case 'qlearning'
        clock_options.agent = 'qlearning';
       [~,~,~,~,mov] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(1:2),params);
    case 'sarsa'
        clock_options.agent = 'sarsa';
        [~,~,~,~,mov] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(1:2),params);
end