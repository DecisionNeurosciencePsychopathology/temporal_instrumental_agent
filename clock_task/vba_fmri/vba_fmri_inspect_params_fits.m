% plot SCEPTIC-SM parameters from VBA_MFX
load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/decay/vba_mfx_results_psub.mat');

% load('/Volumes/bek/vba_results/vba_fmri/fixed_vba_mfx_results_psub.mat');

for ct = 1:length(ids)
alpha(ct) = (p_sub{ct,1}.muTheta(1));
if length(p_sub{ct,1}.muTheta) > 1
gamma(ct) = (p_sub{ct,1}.muTheta(2));
end
beta(ct) = (p_sub{ct,1}.muPhi(1));
end

r1 = round(corr(alpha',gamma'),2);
r2 = round(corr(alpha',beta'),2);


clf; figure(1)
subplot(3,2,1)
scatter(alpha,gamma); xlabel('alphaGauss');ylabel('decayGauss'); text(-3,0,['r = ', num2str(r1)])
subplot(3,2,3); 
hist(alpha,40); xlabel('alphaGauss'); 
subplot(3,2,5)
hist(gamma,40); xlabel('decayGauss')
subplot(3,2,2)
scatter((alpha),(beta));xlabel('alpha');ylabel('beta'); text(-2,3,['r = ', num2str(r2)])
subplot(3,2,4)
hist(beta,40); xlabel('beta')
subplot(3,2,6)
hist(sig(gamma),40); xlabel('decayNative')

% BMC

% decay = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/decay/vba_mfx_results_ogroup.mat');
% fixed = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/fixed/vba_mfx_results_ogroup.mat');

%% test the matrix

%sELECTIVE wIDE fACTORIZED

% wrong, needs to be rerun
% swf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_wide_factorized_vba_mfx_results_ogroup.mat');

swu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_wide_unfactorized_vba_mfx_results_ogroup.mat');
% uniform
uwf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_wide_factorized_vba_mfx_results_ogroup.mat');
uwu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_wide_unfactorized_vba_mfx_results_ogroup.mat');
% uNIFORM nARROW fACTORIZED
unf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_narrow_factorized_vba_mfx_results_ogroup.mat');
unu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uniform_narrow_unfactorize_vba_mfx_L.mat');

% test_swf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_wide_factorized_vba_mfx_results_osub.mat');
% test_swu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_wide_unfactorized_vba_mfx_results_osub.mat');
% 
% test_unf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_narrow_factorized_vba_mfx_results_osub.mat');
% 
% test_snu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_selective_narrow_unfactorize_vba_mfx_results_osub.mat');
% test_snu.o_sub{1}.options.inF
test_unf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_narrow_factorized_vba_mfx_results_osub.mat');
test_unf.o_sub{1}.options.inF

test_uwu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_wide_unfactorized_vba_mfx_results_osub.mat');
test_uwu.o_sub{1}.options.inF
%% snu is in fact snf

%% need to run
swf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_selective_wide_factorize_vba_mfx_L.mat');
%%

snf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_narrow_factorized_vba_mfx_L.mat');
snu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_selective_narrow_unfactorize_vba_mfx_L.mat');
uwu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uniform_wide_unfactorize_vba_mfx_L.mat');
Lswf = swf.F;
Lswu = swu.o_group.within_fit.F;
Lsnu  = snu.F;
Lsnf = snf.F;
Luwf = uwf.o_group.within_fit.F;
Luwu = uwu.F;
% Luwu = uwu.o_group.within_fit.F;
Lunf = unf.o_group.within_fit.F;
Lunu = unu.F;

% % decay = load('/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/decay/vba_mfx_results_ogroup.mat');
% fixed = load('/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/fixed/vba_mfx_results_ogroup.mat');
% uniform = load('/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/mmclock_fmri_decay_uniform_vba_mfx_results_ogroup.mat');
% wide.o_group = o_group;
% % % in this case uniform was just run
% % uni.o_group = o_group;
% 
% % Ldecay = decay.o_group.within_fit.F;
% % outmat is an object form Michael
% 
% % Ldecay_unfactorized = outmat(:,1)';
% Lfixed = fixed.o_group.within_fit.F;
% Luni = uniform.o_group.within_fit.F;
% Lwide = wide.o_group.within_fit.F;
clear options;
options.modelNames = {'SWF';'SWU';'SNF';'SNU';'UWF';'UWU';'UNF';'UNU'};
% sel/uniform
options.families{1} = 1:4;
options.families{2} = 5:8;

% wide/narrow
% options.families{1} = [1,2,5]; % wide
% options.families{2} = [3:4,6:7]; % narrow
% 
% % factorized
options.families{1} = [1,3,5,7]; % factorized
options.families{2} = [2,4,6,8]; % not


[p,o] = VBA_groupBMC([Lswf; Lswu; Lsnf; Lsnu; Luwf;Luwu;Lunf;Lunu],options);

% ignore un-factorized models as theoretically insubstantial
clear options;
options.modelNames = {'SWF';'SNF';'UWF';'UNF'};
% sel/uniform
options.families{1} = 1:2;
options.families{2} = 3:4;
[p,o] = VBA_groupBMC([Lswf; Lsnf; Luwf;Lunf],options);

%
% just the two top contenders -- BOR ns
clear options;
options.modelNames = {'SWF';'UNF'};
[p,o] = VBA_groupBMC([Lswf; Lunf],options);

