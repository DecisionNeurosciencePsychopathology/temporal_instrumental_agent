function run_models(models_to_run)
%This script is designed for single runs of which model you'd like to see
%the input vairbale models_to_run is an array such as [1 2 3 4 5] where 
%1=SKEPTIC
%2=LOGISTIC
%3=GRW
%4=QLEARNING
%5=SARSA


%set up conditional variables
load('s.mat');
load('mIEV.mat')
load('mDEV.mat')
load('mQUADUP.mat')
load('clock_options.mat')

%Set seeds and params here
runseed = 55;
rng(runseed);
seeds = [98 83 66 10];
ntrials=200;
nbasis = 24;
ntimesteps = 500;

%Set diagnosis plots here
trial_plots = 0;
clock_options.diag = trial_plots;

%Set contingency here
m = mIEV; %Set reward matrix can be mIEV mDEV or mQUADUP

%Set if there is a reversal here
rev = 0; %reversal boolean, however m must be either mIEV to mDEV
clock_options.reversal_go=rev;

%Set permuted matrix here
vperm_run = randperm(length(mIEV.sample));
m.vperm_run = vperm_run;


%All models are using optimal params from fitting all conditions via
%Genetic algorithm
for val = 1:length(models_to_run)
    run_mdl = models_to_run(val);
    if run_mdl==1
        %Just run Skeptic
        [cost_SKEPTIC]=clock_logistic_operator_kalman_optimize(s.kalmanSKEPTIC.opt_params(4,:), seeds, m.name, ntrials, nbasis, ntimesteps, trial_plots, 1, m, rev);
    elseif run_mdl==2
        %Just run Logistic
        [cost_Logistic]=clock_logistic_operator_kalman_optimize([s.kalmanLogistic.opt_params(4,:) 0 .2], seeds, m.name, ntrials, nbasis, ntimesteps, trial_plots, 0, m, rev);
    elseif run_mdl==3
        %Just run GRW
        [cost_GRW]=clock_logistic_operator_kalman_optimize([s.kalmanGRW.opt_params(4,:) .9 .2], seeds, m.name, ntrials, nbasis, ntimesteps, trial_plots, 0, m, rev);
    elseif run_mdl==4
        %Just run Q-Learning
        [cost_Q] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(1:2),s.qlearning.opt_params(4,:));
    elseif run_mdl==5
        %Just run SARSA
        [cost_SARSA] = ClockWalking_3D_discountedEv_optimize(clock_options,m, seeds(1:2),s.sarsa.opt_params(4,:));
    end
end
