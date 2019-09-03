function  [fx] = h_sceptic_kalman(x_t, theta, u, inF)
% evolution function of 'kalman' flavors of SCEPTIC
% IN:
%   - x_t : basis values/heights (nbasis x 1)
%   - theta : Model dependent --
%             Fixed uv: theta(1) = tau; theta(2) = alpha;
%             kalman_uv_sum: theta(1) = tau;
%             kalman_softmax: NULL;
%             NOTE: if fit propspread should always be the last parameter of theta
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: evolved basis values/heights (nbasis x 1)

%transform from Gaussian posterior to relevant parameter distributions
theta = transform_theta(theta, inF);

alpha = theta(1); %learning rate (transformed from Gaussian to 0..1)

hidden_state_index=1:inF.hidden_states*inF.nbasis; %total number of hidden states (inF.hidden_states is the number of state vectors)
hidden_state_index = reshape(hidden_state_index, inF.nbasis, inF.hidden_states); %3 x nbasis here

%Define hidden states
%Value should always be the first while uncertainey should always be the second.
w=x_t(hidden_state_index(:,1)); %Value
sigma = x_t(hidden_state_index(:,2)); %posterior uncertainty at each basis

ntimesteps=inF.ntimesteps;
v=w*ones(1,ntimesteps) .* inF.gaussmat;
v_func = sum(v);

% alpha = 1./(1+exp(-theta(1)));
if inF.fit_propspread
  prop_spread = 1./(1+exp(-theta(theta_idx))); %0..1 SD of Gaussian eligibility as proportion of interval
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

%if we are at a run boundary in a u-resetting model, then re-initialize sigma at sigma_noise
if length(u) > 2 && u(3) == 1 && inF.u_run_reset == 1
  sigma(:) = sigma_noise;
end

%compute gaussian spread function with w = rts(i) and sigma based on free param prop_spread
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
e = sum(repmat(elig,inF.nbasis,1).*inF.gaussmat_trunc, 2);

%1) compute prediction error, scaled by eligibility trace
if inF.function_wise_pe
  %%clunky to put this inside observation function
  rtrnd=round(rt);
  if rtrnd < 1
    rtrnd = 1;
  elseif rtrnd > ntimesteps
    rtrnd = ntimesteps;
  end
  delta = e*(reward - v_func(rtrnd));
else
  delta = e.*(reward - w);
end

%Compute the Kalman gains for the current trial (potentially adding process noise)
k = (sigma)./(sigma + sigma_noise);

%Update posterior variances on the basis of Kalman gains
sigma_new = (1 - e.*k).*(sigma);

if strcmpi(inF.model, 'fixed_uv') %this is the u aversion variant
  w_new = w + alpha.*delta; %update values based on fixed learning rate
elseif strcmpi(inF.model, 'fixed_uv_baked')
  %this is the u aversion variant that bakes in uncertainty
  %following the codebase in vba/, I am using the evolved sigmas
  
  tau = theta(2);   %uncertainty approach/avoidance
  w = w + alpha.*delta; %update values based on fixed learning rate
  w_new = w + tau.*sigma_new; %mix together value and uncertainty according to tau *in the representation*
else
  %Standand kalman evolution using dynamic gain as learning rate
  w_new = w + k.*delta;
end

fx = [w_new; sigma_new; delta]; %evolved value, uncertainty, and PE estimates