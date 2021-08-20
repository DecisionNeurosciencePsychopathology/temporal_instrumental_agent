function  [fx] = h_sceptic_exp_compress_fmri(x_t, theta, u, inF)
% evolution function of SCEPTIC model with fixed learning rate and exponential compression of values
%
% The basis weights are the hidden states that evolve with learning.
%
% IN:
%   - x_t : basis values/heights (nbasis x 1)
%   - theta : theta(1) = prop_spread; theta(2) = alpha;
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: evolved basis values/heights (nbasis x 1)

%transform from Gaussian posterior to relevant parameter distributions
theta = transform_theta(theta, inF);

alpha = theta(1); %learning rate (transformed from Gaussian to 0..1)
phi = theta(2); %exponential compression

hidden_state_index=1:inF.hidden_states*inF.nbasis; %total number of hidden states (inF.hidden_states is the number of state vectors)
hidden_state_index = reshape(hidden_state_index, inF.nbasis, inF.hidden_states); %3 x nbasis here

%Define hidden state vectors: mu, sigma, z
w=x_t(hidden_state_index(:,1)); %Value weight vector
%pe = x_t(hidden_state_index(:,2)); %prediction error (not updated per se -- just tracked)

ntimesteps=inF.ntimesteps;
v=w*ones(1,ntimesteps) .* inF.gaussmat;
v_func = sum(v);

if inF.fit_propspread
  prop_spread = inF.max_prop_spread ./ theta(3); %0..max SD of Gaussian eligibility as proportion of interval
  sig_spread=prop_spread*range(inF.tvec); %determine SD of spread function in time units (not proportion)
  
  %if prop_spread is free, then refspread must be recomputed to get AUC of eligibility correct
  refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));
  
  %temporarily re-break the code to check replicability
  %sig_spread = 1 ./ theta(3);
  %refspread = inF.refspread;
else
  sig_spread = inF.sig_spread; %precomputed sig_spread based on fixed prop_spread (see setup_rbf.m)
  refspread = inF.refspread; %precomputed refspread based on fixed prop_spread
end

rt = u(1); %response time 0..40, but can be fractional
reward = u(2); %obtained reward

%compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
elig = gaussmf(inF.tvec, [sig_spread, rt]);

% make up-sampled to 100 Hz space for dirac delta
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

% introduce compression
if inF.exp_variant == 1
  %calculate PE by eligibility-scaled reward - expectation
  delta = e.*(reward - w);

  %update weights by PE
  w_new = w + alpha.*delta;

  % introduce exponential compression
  w_new = exp_renorm(w_new, phi);
elseif inF.exp_variant == 2
  % compress basis weights
  w_compress = exp_renorm(w, phi);

  % calculate PE against compressed weights
  delta = e.*(reward - w_compress);

  % update weights as compressed values plus PE against compressed values
  w_new = w_compress + alpha.*delta;
end

fx=[w_new; delta; w];

end

%small worker function to exponentiate and renormalize onto the same scale
function wtrans = exp_renorm(x, phi)
  if nargin < 1, phi=1; end
  %if x is essentially constant, do not transform
  if all(abs(x - x(1)) < .001)
    %if std(x) < .001 %computationally slower
    wtrans=x;
  else
    trans = exp((x-max(x))/phi); %avoid floating point overflow
    wtrans = (max(x) - min(x))*(trans - min(trans))/(max(trans) - min(trans)) + min(x);
  end
end