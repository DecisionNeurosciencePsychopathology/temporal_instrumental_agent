clock_options_generator('IEV')
load('clock_options.mat')
options.decayflag=0;
options.epsilon=.08;
ClockWalking(options)
