%  runs group model comparison on SCEPTIC log model evidence data from sceptic_fit_group_vba.m

% clear variables;
close all;

cd('/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc/');

% cd('/Users/localadmin/code/clock_smoothoperator/clock_task/vba/results');
load('group_model_comparison_L');
% cd('/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc/');

sceptic.L = L;
sceptic.modelnames = modelnames;
%% Define 'fixed' as the first family, 'kalman', as the second
options.families{1} = find(strcmpi(sceptic.modelnames, 'fixed'));
options.families{2} = find(~strcmpi(sceptic.modelnames, 'fixed'));

[posterior,out] = VBA_groupBMC(sceptic.L,options);

%% save output figure
h = figure(1);
savefig(h,sprintf('VBA_group_BMC_nbasis16_nsteps50'))
