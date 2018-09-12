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

decay = load('/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/decay/vba_mfx_results_ogroup.mat');
fixed = load('/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/fixed/vba_mfx_results_ogroup.mat');
uniform = load('/Volumes/bek/Box Sync/skinner/projects_analyses/SCEPTIC/mfx_analyses/decay/vba_mfx_results_ogroup.mat');

% in this case uniform was just run
uni.o_group = o_group;

Ldecay = decay.o_group.within_fit.F;
% outmat is an object form Michael

% Ldecay_unfactorized = outmat(:,1)';
Lfixed = fixed.o_group.within_fit.F;
Luni = uni.o_group.within_fit.F;

[p,o] = VBA_groupBMC([Lfixed;Ldecay;Luni]);