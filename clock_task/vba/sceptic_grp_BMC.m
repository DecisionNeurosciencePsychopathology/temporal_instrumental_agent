%  runs group model comparison on SCEPTIC log model evidence data from sceptic_fit_group_vba.m

% clear variables;
close all;

cd('/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc/');

% cd('/Users/localadmin/code/clock_smoothoperator/clock_task/vba/results');
load('group_model_comparison_L');
% cd('/Users/localadmin/Google Drive/skinner/SCEPTIC/subject_fitting/vba_results/group_bmc/');

%% really manual 2-model comparison
L(1,:) = L_tight_tau;
L(2,:) = L_kalman_uv_sum;
models = {'tight tau' 'no aversion'};

sceptic.L = L;
sceptic.modelnames = models;


% %% Define 'fixed' as the first family, 'kalman', as the second
% options.families{1} = find(strcmpi(sceptic.modelnames, 'fixed'));
% options.families{2} = find(~strcmpi(sceptic.modelnames, 'fixed'));
% 
% [posterior,out] = VBA_groupBMC(sceptic.L,options);

%% Define 'fixed' as the first family, 'kalman', as the second

[posterior,out] = VBA_groupBMC(sceptic.L);


%% save output figure
h = figure(4);
savefig(h,sprintf('VBA_group_BMC_nbasis4_nsteps10_m1=aversion_supertight_tau_m2=aversion_tight_tau_m3_no_aversion'))
