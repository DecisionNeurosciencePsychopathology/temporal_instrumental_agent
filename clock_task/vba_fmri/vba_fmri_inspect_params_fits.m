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
swf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_wide_factorized_vba_mfx_results_ogroup.mat');
swu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/selective/mmclock_fmri_decay_wide_unfactorized_vba_mfx_results_ogroup.mat');
% uniform
uwf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_wide_factorized_vba_mfx_results_ogroup.mat');
uwu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_wide_unfactorized_vba_mfx_results_ogroup.mat');
% uNIFORM nARROW fACTORIZED
unf = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uni_narrow_factorized_vba_mfx_results_ogroup.mat');
unu = load('~/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/uniform/mmclock_fmri_decay_uniform_narrow_unfactorized_vba_mfx_results_ogroup.mat');

Lswf = swf.o_group.within_fit.F;
Lswu = swu.o_group.within_fit.F;
Luwf = uwf.o_group.within_fit.F;
Luwu = uwu.o_group.within_fit.F;
Lunf = unf.o_group.within_fit.F;
Lunu = unu.o_group.within_fit.F;

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

[p,o] = VBA_groupBMC([Lswf;Lswu;Luwf;Luwu;Lunf;Lunu]);