%  runs group model comparison on SCEPTIC log model evidence data from sceptic_fit_group_vba.m

clear variables;
close all;

cd('/Users/localadmin/code/clock_smoothoperator/clock_task/vba/results');
load('SCEPTIC_6models_011916');

%% Define 'fixed' as the first family, 'kalman', as the second
options.families{1} = find(strcmpi(sceptic.modelnames, 'fixed'));
options.families{2} = find(~strcmpi(sceptic.modelnames, 'fixed'));

[posterior,out] = VBA_groupBMC(sceptic.L,options);
