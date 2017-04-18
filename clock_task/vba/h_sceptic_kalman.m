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


nbasis = inF.nbasis; %Grab basis numbers
hidden_state_index=1:inF.hidden_state*nbasis;
hidden_state_index = reshape(hidden_state_index,nbasis,inF.hidden_state);

%Find length of theta. Depending on the autocorrelation and fitting
%prop_spread we'll need to know this.
theta_idx = length(theta);

%Set the proper paramters NOTE gamma is first then phi for sigma volatility
%model.
if inF.kalman.fixed_uv, alpha = 1./(1+exp(-theta(2))); end %Alpha should always be two as tau should always be 1

%% allow uncertainty aversion in UV_sum?
if inF.kalman.kalman_uv_sum ||  inF.kalman.fixed_uv
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

try
    if strcmp(inF.autocorrelation,'choice_tbf');
        choice = x_t(hidden_state_index(:,end)); %Basis fx's for choice
        theta_idx = theta_idx-1; %Just in casse we want to fit prop_spread
    end
catch
end

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
fx = zeros(size(x_t));

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

%%%%Additional calculations (entropy, v_func,...)%%%%
H = calc_shannon_H(( mu / sum( mu ) )); %Entropy
ntimesteps = inF.ntimesteps;
gaussmat=inF.gaussmat;
v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
v_func = sum(v);

%1) compute prediction error, scaled by eligibility trace
if inF.total_pe %Niv version
    rnd_rt = round(rt);
    
    if rnd_rt==0
        rnd_rt=1;
    elseif rnd_rt>40
        rnd_rt=40;
    end
    delta = e*(reward - v_func(round(rnd_rt)));
else
    delta = e.*(reward - x_t(1:nbasis));
end

%Compute the Kalman gains for the current trial (potentially adding process noise)
k = (sigma)./(sigma + sigma_noise);

%Update posterior variances on the basis of Kalman gains, what used to be
%SIGMA
sigma = (1 - e.*k).*(sigma);
fx(hidden_state_index(:,2)) = sigma;

%I think this is correct [tau * hidden state value + (1-tau) * hidden state uncertainty]
if inF.kalman.kalman_uv_sum
    %Update the value
    mu = mu + k.*delta;
    if inF.u_aversion
        fx(hidden_state_index(:,1))= mu + tau.*sigma; %mix together value and uncertainty according to tau
    else
        fx(hidden_state_index(:,1))=tau.*mu + (1-tau).*sigma; %mix together value and uncertainty according to tau
    end
    
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


% The try catches are in place for refitting purposes of older datasets
% that never had the following options

%track pe as a hidden state
try
    if inF.track_pe
        fx(hidden_state_index(:,end))=delta;
    end
catch
end


%Entropy update to hidden states
try
    if inF.entropy
        H = is_nan_or_inf(H);
        %max_value = is_nan_or_inf(max_value);
        max_value = find_max_value_in_time(v_func);
        fx(hidden_state_index(end)+1) = H;
        fx(hidden_state_index(end)+2) = max_value;
    end
catch
end

try
    if strcmp(inF.autocorrelation,'choice_tbf')
        choice_decay = 1./(1+exp(-theta(end)));
        fx(hidden_state_index(:,end)) = choice_decay.*choice + e;
    end
catch
end

%FOR uv sum we need to to rewrite fx(1:nbasis) to be the uv sum analog to
%the old clock sceptic script


function out = is_nan_or_inf(var)
if isnan(var) || isinf(var)
    out=0;
else
    out=var;
end

function time_point = find_max_value_in_time(val)
max_val  = max(val); %Get the max value
time_point = find(max_val==val); %Find the max
if length(time_point)>1 %If there are duplicates randomly select one
    time_point = randsample(time_point,1);
end


