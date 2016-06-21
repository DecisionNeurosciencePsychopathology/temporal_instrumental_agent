function  [fx] = h_sceptic_fixed_decay(x_t, theta, u, inF)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : basis values/heights (nbasis x 1)
%   - theta : theta(1) = prop_spread; theta(2) = alpha;
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: evolved basis values/heights (nbasis x 1)

alpha = 1./(1+exp(-theta(1))); %LR: 0..1
gamma = 1./(1+exp(-theta(2))); %Decay: 0..1
if inF.fit_propspread
    prop_spread = 1./(1+exp(-theta(3))); %0..1 SD of Gaussian eligibility as proportion of interval
    sig_spread=prop_spread*range(inF.tvec); %determine SD of spread function in time units (not proportion)
    
    %if prop_spread is free, then refspread must be recomputed to get AUC of eligibility correct
    refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));
else
    sig_spread = inF.sig_spread; %precomputed sig_spread based on fixed prop_spread (see setup_rbf.m)
    refspread = inF.refspread; %precomputed refspread based on fixed prop_spread
end

rt = u(1);
reward = u(2);
if inF.fit_nbasis
    %% convert normally distributed theta(2) to discrete uniform number of bases
    nbasis_cdf = cdf('Normal',theta(2),inF.muTheta2, inF.SigmaTheta2);
    nbasis = unidinv(nbasis_cdf,inF.maxbasis);
else
    nbasis = inF.nbasis;
end

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
e = sum(repmat(elig,nbasis,1).*inF.gaussmat_trunc, 2);

%1) compute prediction error, scaled by eligibility trace
delta = e.*(reward - x_t);

%% introduce decay
decay = -gamma.*(1-e).*x_t;

fx = x_t + alpha.*delta + decay;


