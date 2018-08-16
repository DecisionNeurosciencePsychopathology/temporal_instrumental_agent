% plot SCEPTIC-SM parameters from VBA_MFX

load('/Volumes/bek/vba_results/vba_fmri/fixed_vba_mfx_results_psub.mat');

for ct = 1:length(id)
group_alpha(ct) = (p_sub{ct,1}.muTheta(1));
if length(p_sub{ct,1}.muTheta(2)) > 1
group_decay(ct) = (p_sub{ct,1}.muTheta(2));
end
end

r = round(corr(group_alpha',group_decay'),2);


clf; figure(1)
subplot(3,2,1)
scatter(group_alpha,group_decay); xlabel('alphaGauss');ylabel('decayGauss'); text(-5,0,['r = ', num2str(r)])
subplot(3,2,3); 
hist(group_alpha,40); xlabel('alphaGauss'); 
subplot(3,2,5)
hist(group_decay,40); xlabel('decayGauss')
subplot(3,2,2)
scatter(sig(group_alpha),sig(group_decay));xlabel('alphaNative');ylabel('decayNative')
subplot(3,2,4)
hist(sig(group_alpha),40); xlabel('alphaNative')
subplot(3,2,6)
hist(sig(group_decay),40); xlabel('decayNative')

