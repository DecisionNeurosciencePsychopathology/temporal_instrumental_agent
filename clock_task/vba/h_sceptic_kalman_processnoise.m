function  [fx] = h_sceptic_kalman_processnoise(x_t, theta, u, inF)
% This is the evolution function of the kalman processnoise version of the
% sceptic model.
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : basis values/heights (nbasis x 1)
%   - theta : theta will equal omega, discrim, tau, gamma, or phi.  theta(end) = prop_spread
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: evolved basis values/heights (nbasis x 1)

%FOR reference
omega=0;                    %if not a PE-enhanced process noise model, remove this influence
gamma=0;                    %zero gamma and phi for all models except kalman_sigmavolatility
phi=0;
tradeoff=0;                 %parameter only applies to kalman_uv_logistic
z=0;

%Set the proper paramters NOTE gamma is first then phi for sigma volatility
%model.
if inF.kalman.processnoise, omega = 1./(1+exp(-theta(1))); end
if inF.kalman.sigmavolatility
    gamma = 1./(1+exp(-theta(1)));
    inF.kalman.phi, phi = 1./(1+exp(-theta(2)));
end
if inF.kalman.uv_logistic, tradeoff = 1./(1+exp(-theta(1))); end




% alpha = 1./(1+exp(-theta(1)));
if inF.fit_propspread
    prop_spread = 1./(1+exp(-theta(end))); %0..1 SD of Gaussian eligibility as proportion of interval
    sig_spread=prop_spread*range(inF.tvec); %determine SD of spread function in time units (not proportion)
    
    %if prop_spread is free, then refspread must be recomputed to get AUC of eligibility correct
    refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));
else
    sig_spread = inF.sig_spread; %precomputed sig_spread based on fixed prop_spread (see setup_rbf.m)
    refspread = inF.refspread; %precomputed refspread based on fixed prop_spread
end

rt = u(1);
reward = u(2);
sigma_noise = inF.sigma_noise;
fx = zeros(size(x_t));
% There is no other argument in kalman softamax but prop_spread
%if inF.fit_nbasis
%% convert normally distributed theta(2) to discrete uniform number of bases
%    nbasis_cdf = cdf('Normal',theta(2),inF.muTheta2, inF.SigmaTheta2);
%    nbasis = unidinv(nbasis_cdf,inF.maxbasis);
%else
nbasis = inF.nbasis;
%end
% ntimesteps = inF.ntimesteps;

%refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));

%compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
% elig = gaussmf(inF.tvec, [inF.sig_spread, rt]);
elig = gaussmf(inF.tvec, [sig_spread, rt]);


%compute sum of area under the curve of the gaussian function
auc=sum(elig);

%divide gaussian update function by its sum so that AUC=1.0, then rescale to have AUC of a non-truncated basis
%this ensures that eligibility is 0-1.0 for non-truncated update function, and can exceed 1.0 at the edge.
%note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
%will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
elig=elig/auc*refspread;

%compute the intersection of the Gaussian spread function with the truncated Gaussian basis.
%this is essentially summing the area under the curve of each truncated RBF weighted by the truncated
%Gaussian spread function.
e = sum(repmat(elig,nbasis,1).*inF.gaussmat_trunc, 2);

%1) compute prediction error, scaled by eligibility trace
delta = e.*(reward - x_t(1:nbasis));



%Kalman variants of learning and uncertainty

%Changing Kalman variance a posteriori uses the elig*gain approach: [1 - k(ij)*elig(ij)]*sigma(ij)
%this would only allow a 1.0 update*kalman gain for basis functions solidly in the window and a decay in diminishing
%variance as the basis deviates from the timing of the obtained reward.


% sigma = sigma_noise; BUT SET SIGMA_NOISE PRIORS FOR EACH SESSION
sigma = x_t(nbasis+1:end);
%THIS will be nothing since there is no omega term...
%MH 8Sep2015: At the moment, we assume zero process noise in the estimated posterior error covariances, sigma_ij.
%To model dynamic/changing systems, try dynamically enhance learning rates by scaling process noise by PE.
Q = omega.*abs(delta); %use abs of PE so that any large surprise enhances effective gain.

%Compute the Kalman gains for the current trial (potentially adding process noise)
k = (sigma + Q)./(sigma + Q + sigma_noise);

%Update posterior variances on the basis of Kalman gains, what used to be
%SIGMA
fx(nbasis+1:end) = (1 - e.*k).*(sigma + z);

%Update reward expectation. AD: Would it be better for the delta to be the difference between the reward
%and the point value estimate at the RT(i)?

%This might be the end for this guy!...
% mu_ij = mu_ij(i,:) + k_ij(i,:).*delta_ij(i,:);

%For this case these drop out as well
%Track smooth estimate of volatility according to unsigned PE history
z = gamma.*z + phi.*abs(sum(delta));

%Uncertainty is a function of Kalman uncertainties.
% u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat;
% u_func = sum(u_jt); %vector of uncertainties by timestep
% u_it(i+1,:) = u_func;

% What our final for kalman will be is x_t + k_ij * delta so we need to
% compute k
fx(1:nbasis) = x_t(1:nbasis) + k.*delta;






