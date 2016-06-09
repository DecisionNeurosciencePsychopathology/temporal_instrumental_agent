function  [fx] = h_sceptic_kalman_monster(x_t, theta, u, inF)
% This is the evolution function of the kalman monster sceptic model (reversion to prior, volatility, precision weighting)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx] = h_sceptic_kalman_monster(x_t,theta,u,inF)
% IN:
%   - x_t : basis values/heights (nbasis x 1)
%   - theta : theta will equal omega, discrim, tau, gamma, or phi.  theta(end) = prop_spread
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: evolved basis values/heights (nbasis x 1)

%FOR reference
omega=0;                    %if not a PE-enhanced process noise model, remove this influence

nbasis = inF.nbasis; %Grab basis numbers
hidden_state_index=1:inF.hidden_state*nbasis; %total number of hidden states (inF.hidden_state is the number of state vectors)
hidden_state_index = reshape(hidden_state_index,nbasis,inF.hidden_state); %3 x nbasis here

%Define hidden state vectors: mu, sigma, z
mu=x_t(hidden_state_index(:,1)); %Value state vector
sigma = x_t(hidden_state_index(:,2)); %Uncertainty state vector
z = x_t(hidden_state_index(:,3)); %Volatility state vector

%extract thata parameters: kappa, gamma, phi, alpha
kappa = theta(1); %should range from about .1 - 10
gamma = 1./(1+exp(-theta(2))); %AR term: sigmoid transform 0..1
%phi = 1./(1+exp(-theta(3)));
phi = theta(3); %should be Gaussian distributed around 0, small variance (negative would a weird result [PE reduces volatility], but not impossible)
alpha = 1./(1+exp(-theta(4))); %sigmoid for power scaling of uncertainty

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
prior_sigma = inF.priors.muX0(hidden_state_index(:,2)); %prior values for sigma state vector

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
%this is essentially summing the area under the curve of each truncated RBF weighted by the truncated Gaussian spread function.
e = sum(repmat(elig,nbasis,1).*inF.gaussmat_trunc, 2);

%compute prediction error, scaled by eligibility trace
delta = e.*(reward - mu);

%MH 8Sep2015: At the moment, we assume zero process noise in the estimated posterior error covariances, sigma_ij.
%To model dynamic/changing systems, try dynamically enhance learning rates by scaling process noise by PE.
Q = omega.*abs(delta); %use abs of PE so that any large surprise enhances effective gain.

%Compute the Kalman gains for the current trial (potentially adding process noise)
k = (sigma + Q)./(sigma + Q + sigma_noise);

%Update means (mu)
fx(hidden_state_index(:,1)) = mu + k.*delta;

%compute Bayesian forgetting operator
psi = 1./(1 + exp(-kappa.*(z + ((1 - e).*prior_sigma).^alpha)));

%Update posterior variances (sigma) on the basis of reduction in uncertainty versus forgetting and volatility
fx(hidden_state_index(:,2)) = ((1 - e.*k).*sigma).^psi .* prior_sigma.^(1-psi);

%Track smooth estimate of volatility according to unsigned PE history locally instead of globally
fx(hidden_state_index(:,3)) = (1-e).*z + e.*(gamma*z + (prior_sigma./sigma).*phi.*abs(sum(delta)));
%fprintf('fx: %s\n', num2str(fx));
%disp(fx)
if any(isnan(fx)), disp(fx); end
if any(isinf(fx)), disp(fx); end