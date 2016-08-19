clear all
%cd '/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing'

%load('output/optimize_output_sinusoid_fixedLR_softmax.mat');
load('sinusoid_optmat.mat');
%p1=pars{1};
p1=[0.1332, 0.3992, 0.0696];

a_all=initialize_agents_struct;
fixedlr=a_all(2);
fixedlr.nbasis = 24;
fixedlr.ntimesteps = 500;
fixedlr.ntrials = 100;
fixedlr.runseed = 888; %seeds rngs inside runs?

reps=10
totcosts=zeros(1,reps);
for r = 1:reps
    [totcosts(r), perruncosts, seeds] = multirun_clock_sceptic(p1, fixedlr, optmat{1});
end

round(totcosts)
round(costs(1))

fixedlr.init_params=p1;
[Xga,Fga] = GAoptimize(fixedlr, optmat{1});



qlearning = a_all(14);
qlearning.nbasis = 24;
qlearning.ntimesteps = 500;
qlearning.ntrials = 100;
qlearning.runseed = 888; %seeds rngs inside runs?
qlearning.clock_options.episodeCount = qlearning.ntrials;
qlearning.clock_options.ntimesteps = qlearning.ntimesteps/10;

multirun_cliff_walker_optimize(qlearning.init_params, qlearning, optmat{1})