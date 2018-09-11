load vba_mfx_results_pgroup.mat
load vba_mfx_results_psub.mat
load vba_mfx_results_ogroup.mat
load vba_mfx_results_osub.mat

%add VBA to MATLAB's universe
addpath(genpath('/gpfs/group/mnh5174/default/lab_resources/VBA-toolbox'))

o_group.options.DisplayWin=1; %need to set this post hoc for the display to work
VBA_displayMFX(p_sub, o_sub, p_group, o_group);

ns=length(p_sub);

%get log evidence vector
outmat = NaN(ns, 3);
for (nn=1:ns)
  outmat(nn,:) = [o_group.within_fit.F(nn) o_group.within_fit.R2(nn) o_group.within_fit.LLH0(nn)];
end

save('mfx_decay_fixed_x0.mat', 'outmat');

n_theta=2;
n_phi=1;

Theta=NaN(ns, n_theta, 2); %first element is MFX approach, second is typical within subject
Phi=NaN(ns, n_phi, 2); 

for i=1:ns
  %Theta(i,:,1) = sig(p_sub{i}.muTheta); %inv logit transform on alpha and gamma
  %Phi(i,:,1) = exp(p_sub{i}.muPhi); %exponential transform on temperature
  %Theta(i,:,2) = sig(o_group.initVBA.p_sub{i}.muTheta);
  %Phi(i,:,2) = exp(o_group.initVBA.p_sub{i}.muPhi);  
  
  Theta(i,:,1) = p_sub{i}.muTheta; %inv logit transform on alpha and gamma
  Phi(i,:,1) = p_sub{i}.muPhi; %exponential transform on temperature
  Theta(i,:,2) = o_group.initVBA.p_sub{i}.muTheta;
  Phi(i,:,2) = o_group.initVBA.p_sub{i}.muPhi;

end

%alpha: learning rate
figure(2);
subplot(3,2,1);
hist(Theta(:,1,1));
title('alpha: mfx estimate')

subplot(3,2,2);
hist(Theta(:,1,2));
title('alpha: single subject estimate')

%gamma: selective maintenance
subplot(3,2,3);
hist(Theta(:,2,1));
title('gamma: mfx estimate')

subplot(3,2,4);
hist(Theta(:,2,2));
title('gamma: single subject estimate')

%beta: temperature
subplot(3,2,5);
hist(Phi(:,1,1));
title('beta: mfx estimate')

subplot(3,2,6);
hist(Phi(:,1,2));
title('beta: single subject estimate')

figure(3);
scatter(Theta(:,1,2), Theta(:,2,2))
corr(Theta(:,1,2), Theta(:,2,2))

corr(Theta(:,1,1), Theta(:,1,2)) %alpha between mfx and single
corr(Theta(:,2,1), Theta(:,2,2)) %gamma between mfx and single






