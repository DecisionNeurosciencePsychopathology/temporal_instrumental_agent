function  [fx] = h_sceptic_fixed_fmri(x_t, theta, u, inF)
% evolution function of SCEPTIC model with fixed learning rate and decay
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

% factorize = 0;

alpha = 1./(1+exp(-theta(1))); %learning rate (transformed from Gaussian to 0..1)
% gamma = 1./(1+exp(-theta(2))); %decay rate (transformed from Gaussian to 0..1)

hidden_state_index=1:inF.hidden_state*inF.nbasis; %total number of hidden states (inF.hidden_state is the number of state vectors)
hidden_state_index = reshape(hidden_state_index,inF.nbasis,inF.hidden_state); %3 x nbasis here

%Define hidden state vectors: mu, sigma, z
w=x_t(hidden_state_index(:,1)); %Value weight vector
%pe = x_t(hidden_state_index(:,2)); %prediction error (not updated per se -- just tracked)
%decay = x_t(hidden_state_index(:,3)); %decay vector (not updated, just tracked)

ntimesteps=inF.ntimesteps;
v=w*ones(1,ntimesteps) .* inF.gaussmat;
v_func = sum(v);

if inF.fit_propspread
    prop_spread = inF.max_prop_spread ./ (1+exp(-theta(3))); %0..max SD of Gaussian eligibility as proportion of interval
    sig_spread=prop_spread*range(inF.tvec); %determine SD of spread function in time units (not proportion)
    
    %if prop_spread is free, then refspread must be recomputed to get AUC of eligibility correct
    refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));
    
    %temporarily re-break the code to check replicability
    %sig_spread = 1./(1+exp(-theta(3)));
    %refspread = inF.refspread;
else 
    sig_spread = inF.sig_spread; %precomputed sig_spread based on fixed prop_spread (see setup_rbf.m)
    refspread = inF.refspread; %precomputed refspread based on fixed prop_spread
end

rt = u(1); %response time 0..40, but can be fractional
reward = u(2); %obtained reward

%compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
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

%compute prediction error, scaled by eligibility trace
nivpe=0; %using basis-wise V for now
if nivpe
  %%clunky to put this inside observation function
  rtrnd=round(rt);
  if rtrnd < 1
    rtrnd = 1;
  elseif rtrnd > ntimesteps
    rtrnd = ntimesteps;
  end
    
  delta = e.*(reward - v_func(rtrnd));
else
  delta = e.*(reward - w);
end

% % introduce decay
% if factorize
%     decay = -gamma.*alpha.*(1-e).*w;
% else
%     decay = -gamma.*(1-e).*w;
% end
w_new = w + alpha.*delta;

fx=[w_new; delta];

end