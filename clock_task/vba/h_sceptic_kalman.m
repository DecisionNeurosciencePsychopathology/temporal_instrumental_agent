function  [fx] = h_sceptic_kalman(x_t, theta, u, inF)
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

nbasis = inF.nbasis; %Grab basis numbers
hidden_state_index=1:inF.hidden_state*nbasis;
hidden_state_index = reshape(hidden_state_index,nbasis,inF.hidden_state);

%Set the proper paramters NOTE gamma is first then phi for sigma volatility
%model.
if inF.kalman.processnoise, omega = 1./(1+exp(-theta(1))); end
if inF.kalman.kalman_sigmavolatility || inF.kalman.kalman_sigmavolatility_local || inF.kalman.kalman_sigmavolatility_precision
    if inF.no_gamma
        gamma = 1-phi;
    else
        gamma = 1./(1+exp(-theta(2)));
    end
    phi = 1./(1+exp(-theta(1)));
    z = x_t(hidden_state_index(:,3)); %Volatility
end
if inF.kalman.kalman_uv_sum_sig_vol
    gamma = 1./(1+exp(-theta(2)));
    phi = 1./(1+exp(-theta(3)));
    %     if inF.u_aversion
    %         tau = theta(3)./1000; % scale it down a bit
    %     else
    %         tau = 1./(1+exp(-theta(3)-log(inF.sigma_noise)));
    %     end
    %     z = x_t(hidden_state_index(:,3)); %Volatility
end
if inF.kalman.kalman_uv_logistic, tradeoff = 1./(1+exp(-theta(1))); end
if inF.kalman.fixed_uv, alpha = 1./(1+exp(-theta(2))); end %Alpha should always be two as tau should always be 1

%% allow uncertainty aversion in UV_sum?
if inF.kalman.kalman_uv_sum || inF.kalman.kalman_uv_sum_sig_vol || inF.kalman.fixed_uv
    if inF.u_aversion
        tau = theta(1)./1000; % scale it down a bit
    else
        tau = 1./(1+exp(-theta(1)-log(inF.sigma_noise)));
    end
end

%Define hidden states
%Value should always be the first while uncertainey should always be the
%second.
mu=x_t(hidden_state_index(:,1)); %Value
sigma = x_t(hidden_state_index(:,2)); %Uncertainty


% mu=x_t(1:nbasis); %Value
% sigma = x_t(nbasis+1:end); %Uncertainty





% alpha = 1./(1+exp(-theta(1)));
if inF.fit_propspread
    prop_spread = 1./(1+exp(-theta(end)));
else
    prop_spread = inF.sig_spread;
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

%end
% ntimesteps = inF.ntimesteps;

%refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));

%compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
% elig = gaussmf(inF.tvec, [inF.sig_spread, rt]);
elig = gaussmf(inF.tvec, [prop_spread, rt]);


%compute sum of area under the curve of the gaussian function
auc=sum(elig);

%divide gaussian update function by its sum so that AUC=1.0, then rescale to have AUC of a non-truncated basis
%this ensures that eligibility is 0-1.0 for non-truncated update function, and can exceed 1.0 at the edge.
%note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
%will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
elig=elig/auc*inF.refspread;

%compute the intersection of the Gaussian spread function with the truncated Gaussian basis.
%this is essentially summing the area under the curve of each truncated RBF weighted by the truncated
%Gaussian spread function.
e = sum(repmat(elig,nbasis,1).*inF.gaussmat_trunc, 2);




%1) compute prediction error, scaled by eligibility trace
delta = e.*(reward - mu);



%Kalman variants of learning and uncertainty

%Changing Kalman variance a posteriori uses the elig*gain approach: [1 - k(ij)*elig(ij)]*sigma(ij)
%this would only allow a 1.0 update*kalman gain for basis functions solidly in the window and a decay in diminishing
%variance as the basis deviates from the timing of the obtained reward.


%MH 8Sep2015: At the moment, we assume zero process noise in the estimated posterior error covariances, sigma_ij.
%To model dynamic/changing systems, try dynamically enhance learning rates by scaling process noise by PE.
Q = omega.*abs(delta); %use abs of PE so that any large surprise enhances effective gain.

%Compute the Kalman gains for the current trial (potentially adding process noise)
k = (sigma + Q)./(sigma + Q + sigma_noise);

%Update posterior variances on the basis of Kalman gains, what used to be
%SIGMA
sigma = (1 - e.*k).*(sigma + z);
fx(hidden_state_index(:,2)) = sigma;
%
%
% %THIS need to be fx(whatever) to update next timestep
% z = gamma.*z + phi.*abs(sum(delta));


%Uncertainty is a function of Kalman uncertainties.
%This may be un-needed...
%u=sigma'*ones(1,inF.ntimesteps) .* inF.gaussmat;
%u_func = sum(sigma); %vector of uncertainties by timestep



%I think this is correct [tau * hidden state value + (1-tau) * hidden state uncertainty]
if inF.kalman.kalman_uv_sum
    %Update the value
    mu = mu + k.*delta;
    if inF.u_aversion
        fx(hidden_state_index(:,1))= mu + tau.*sigma; %mix together value and uncertainty according to tau
    else
        fx(hidden_state_index(:,1))=tau.*mu + (1-tau).*sigma; %mix together value and uncertainty according to tau
        %fx(1:nbasis)=tau.*mu + (1-tau).*fx(nbasis+1:end); %mix together value and uncertainty according to tau
    end
elseif inF.kalman.kalman_sigmavolatility
    %Update value
    fx(hidden_state_index(:,1)) = mu + k.*delta;
    %Track smooth estimate of volatility according to unsigned PE history
    fx(hidden_state_index(:,3)) = gamma.*z + phi.*abs(sum(delta));
elseif inF.kalman.kalman_sigmavolatility_local
    %Update value
    fx(hidden_state_index(:,1)) = mu + k.*delta;
    
    %Overwrite the sigma term
    fx(hidden_state_index(:,2)) = (1 - e.*k).*sigma + z;
    
    %Track smooth estimate of volatility according to unsigned PE history
    %locally instead of globally
    fx(hidden_state_index(:,3)) = (1-e).*z + e.*(gamma*z + phi.*abs(sum(delta)));
elseif inF.kalman.kalman_sigmavolatility_precision
    %Grab priors
    prior_sigma = inF.priors.muX0(hidden_state_index(:,2));
    
    %Update value
    fx(hidden_state_index(:,1)) = mu + k.*delta;
    
    %Overwrite the sigma term
    fx(hidden_state_index(:,2)) = (1 - e.*k).*sigma + z;
    
    %Track smooth estimate of volatility according to unsigned PE history
    %locally instead of globally
    fx(hidden_state_index(:,3)) = (1-e).*z + e.*(gamma*z + (prior_sigma./sigma)*phi.*abs(sum(delta)));
elseif inF.kalman.kalman_uv_sum_sig_vol
    %Update the value
    mu = mu + k.*delta;
    if inF.u_aversion
        fx(hidden_state_index(:,1))= mu + tau.*sigma; %mix together value and uncertainty according to tau
    else
        fx(hidden_state_index(:,1))=tau.*mu + (1-tau).*sigma; %mix together value and uncertainty according to tau
    end
    fx(hidden_state_index(:,3)) = gamma.*z + phi.*abs(sum(delta));
    
elseif inF.kalman.fixed_uv
    %Track sigma use fixed value update, so we'll need alpha in the mix
    mu = mu + alpha.*delta;
    if inF.u_aversion
        fx(hidden_state_index(:,1))= mu + tau.*sigma; %mix together value and uncertainty according to tau
    else
        fx(hidden_state_index(:,1))=tau.*mu + (1-tau).*sigma; %mix together value and uncertainty according to tau
    end
    
else
    % What our final for kalman will be is x_t + k_ij * delta so we need to
    % compute k
    fx(hidden_state_index(:,1)) = mu + k.*delta;
end

%FOR uv sum we need to to rewrite fx(1:nbasis) to be the uv sum analog to
%the old clock sceptic script




