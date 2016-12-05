function [p, agent] = get_sceptic_parameters(params, agent)
%this function loads parameters from the unnamed parameter vector params into the named structure p.
%this is intended to allow for more flexibility of parameterization variants and to support more effectively the
%use of a fixed beta for any agent using the softmax choice rule. Optimizers pass in an unnamed parameter vector,
%and the elements of this vector depend on the model variant, detracting from simplicity in parameter initialization.

p.omega=0;                    %if not a PE-enhanced process noise model, remove this influence
p.gamma=0;                    %zero gamma and phi for all models except kalman_sigmavolatility
p.phi=0;
p.tradeoff=0;                 %parameter only applies to kalman_uv_logistic

%can pass in agent structure for easier named lookup of parameters
useastruct=0;
fixbeta=0;
if isstruct(agent)
    useastruct=1;
    astruct=agent;
    agent=astruct.name;
    if isfield(astruct, 'fixbeta'), fixbeta=astruct.fixbeta; end %use fixed beta from agent if specified
    if isfield(astruct, 'fixps') && astruct.fixps > 0
        p.prop_spread=astruct.fixps;
    else
        p.prop_spread = params(strcmpi(astruct.parnames, 'prop_spread'));
    end
else
    %if not using agent structure, always assume prop spread is first parameter
    p.prop_spread = params(1);      %proportion of discrete interval over which to spread reward (SD of Gaussian) (0..1)
end


if fixbeta > 0
    %populate fixed beta if specified
    p.beta=fixbeta;
elseif useastruct && any(ismember(astruct.parnames, 'beta'))
    %populate beta parameter if agent uses it
    p.beta = params(strcmpi(astruct.parnames, 'beta'));
end

%populate free parameters for the requested agent/model
%note: Kalman filter does not have a free learning rate parameter.
if strcmpi(agent, 'fixedLR_softmax')
    %Learning: Bush-Mosteller fixed learning rate applied to PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: softmax function of value to select among alternatives
    if useastruct
        p.alpha = params(strcmpi(astruct.parnames, 'alpha'));
    else
        p.beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)
        p.alpha = params(3);          %learning rate (0..1)
    end
elseif strcmpi(agent, 'fixedLR_egreedy')
    %Learning: Bush-Mosteller fixed learning rate applied to PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: Epsilon-greedy with hard max of value for exploit and random uniform explore
    if useastruct
        p.epsilon = params(strcmpi(astruct.parnames, 'epsilon'));
        p.alpha = params(strcmpi(astruct.parnames, 'alpha'));
    else
        p.epsilon = params(2);        %proportion of exploratory choices (0..1)
        p.alpha = params(3);          %learning rate (0..1)
    end
    
elseif strcmpi(agent, 'fixedLR_egreedy_grw')
    %Learning: Bush-Mosteller fixed learning rate applied to PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: Epsilon-greedy with hard max of value for exploit and GRW exploration
    if useastruct
        p.epsilon = params(strcmpi(astruct.parnames, 'epsilon'));
        p.alpha = params(strcmpi(astruct.parnames, 'alpha'));
        p.sig_grw = params(strcmpi(astruct.parnames, 'sig_grw'));
    else
        p.epsilon = params(2);        %proportion of exploratory choices (0..1)
        p.alpha = params(3);          %learning rate (0..1)
        p.sig_grw = params(4);        %SD of GRW (0..1 interval proportion) rescaled wrt the observed interval
    end
    
    %determine whether to strategic explore versus GRW (not currently used)
    %rng(exptype_rng_seed);
    %explore_type_rand=rand(ntrials,1);
    %exptype_rng_state=rng; %no need to save if this is not re-used.
elseif strcmpi(agent, 'asymfixedLR_softmax')
    %Learning: Bush-Mosteller separate fixed learning rates for PE+ and PE- according to radial basis; no representation of uncertainty
    %Choice: softmax function of value to select among alternatives    
    if useastruct
        p.alpha = params(strcmpi(astruct.parnames, 'alpha'));
        p.rho = params(strcmpi(astruct.parnames, 'rho'));
    else
        p.beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)
        p.alpha = params(3);          %learning rate for PE+ (0..1)
        p.rho = params(4);            %learning rate for PE- (0..1)
    end
elseif strcmpi(agent, 'kalman_softmax')
    %Learning: kalman filter, gain initialized to 0.5 with Bayesian update (roughly exponential decay); uncertainty is tracked, but not used
    %Choice: softmax function of value to select among alternatives (no uncertainty)
    if ~useastruct
        p.beta = params(2);           %temperature parameter scaling preference among alternatives in softmax (.001..2)
    end
elseif strcmpi(agent, 'kalman_processnoise')
    %Learning: kalman filter, gain initialized to 0.5 with Bayesian update (roughly exponential decay); uncertainty is tracked, but not used
    %Choice: softmax function of value to select among alternatives
    %Additional: PEs boost process noise (Q), effectively enhancing learning rate when surprising events occur.
    if useastruct
        p.omega = params(strcmpi(astruct.parnames, 'omega'));
    else
        p.beta  = params(2);          %temperature parameter scaling preference among alternatives in softmax (.001..2)
        p.omega = params(3);          %scaling of process noise Q by PE.
    end
elseif ismember(agent, {'kalman_sigmavolatility', 'kalman_sigmavolatility_local', 'kalman_sigmavolatility_local_precision'})
    %Learning: kalman filter, gain initialized to 0.5 with Bayesian update (roughly exponential decay); no representation of uncertainty
    %Choice: softmax function of value to select among alternatives
    %Additional: Prediction errors increase uncertainty (sigma) by a smooth function of PEs according to a second learning rate, phi.
    if useastruct
        p.phi   = params(strcmpi(astruct.parnames, 'phi'));
        p.gamma = params(strcmpi(astruct.parnames, 'gamma'));
    else
        p.beta  = params(2);          %temperature parameter scaling preference among alternatives in softmax (.001..2)
        p.phi   = params(3);          %scales additional noise added to sigma update according to smooth function of PEs.
        p.gamma = params(4);          %decay factor for volatility (0..1)
    end
elseif strcmpi(agent, 'kalman_uv_logistic')
    %old logistic explore/exploit choice rule
    if useastruct
        p.tradeoff = params(strcmpi(astruct.parnames, 'tradeoff'));
        p.discrim  = params(strcmpi(astruct.parnames, 'discrim'));
    else
        p.tradeoff = params(2); %point at which agent is indifferent between explore and exploit
        p.discrim  = params(3); %discrimination/slope of logistic function for explore/exploit sampling (0.01..100)
    end
elseif strcmpi(agent, 'kalman_uv_sum')
    if useastruct
        p.tau  = params(strcmpi(astruct.parnames, 'tau'));
    else
        p.beta = params(2);
        p.tau  = params(3); % tau from Greek ???? -- value, price, cf ???? as trophys in Homer
    end
elseif strcmpi(agent, 'fixedLR_kl_softmax')
    if useastruct
        p.alpha = params(strcmpi(astruct.parnames, 'alpha'));
        p.kappa = params(strcmpi(astruct.parnames, 'kappa'));
        p.lambda = params(strcmpi(astruct.parnames, 'lambda'));
    else
        p.beta = params(2);
        p.alpha = params(3);   %learning rate for value
        p.kappa = params(4);   %PE+ tilt scaling parameter
        p.lambda = params(5);  %PE- tilt scaling parameter
    end
elseif strcmpi(agent, 'kalman_kl_softmax')
    if useastruct
        p.kappa = params(strcmpi(astruct.parnames, 'kappa'));
        p.lambda = params(strcmpi(astruct.parnames, 'lambda'));
    else
        p.beta = params(2);
        p.kappa = params(3);
        p.lambda = params(4);
    end
elseif strcmpi(agent, 'kalman_processnoise_kl')
    if useastruct
        p.omega = params(strcmpi(astruct.parnames, 'omega'));
        p.kappa = params(strcmpi(astruct.parnames, 'kappa'));
        p.lambda = params(strcmpi(astruct.parnames, 'lambda'));
    else
        p.beta = params(2);
        p.omega = params(3);  %scaling of process noise by PE.
        p.kappa = params(4);  %PE+ tilt scaling parameter
        p.lambda = params(5); %PE- tilt scaling parameter
    end
elseif strcmpi(agent, 'kalman_uv_sum_kl')
    if useastruct
        p.tau  = params(strcmpi(astruct.parnames, 'tau'));
        p.kappa = params(strcmpi(astruct.parnames, 'kappa'));
        p.lambda = params(strcmpi(astruct.parnames, 'lambda'));
    else
        p.beta = params(2);
        p.tau = params(3);
        p.kappa = params(4);
        p.lambda = params(5);
    end
elseif ismember(agent, {'kalman_uv_sum_discount', 'kalman_uv_sum_negtau'})
    if useastruct
        p.tau  = params(strcmpi(astruct.parnames, 'tau'));
    else
        p.beta = params(2);
        p.tau = params(3);
        p.kappa = params(4);
        p.lambda = params(5);
    end
elseif ismember(agent, {'fixedLR_decay', 'fixedLR_decay_random'})
    if useastruct
        p.alpha  = params(strcmpi(astruct.parnames, 'alpha'));
        p.gamma = 1/(1+exp(-params(strcmpi(astruct.parnames, 'gamma')))); %exponentiate to convert to 0-1 scaling
    else
        p.beta = params(2);
        p.alpha = params(3);
        p.gamma = params(4);
    end
elseif ismember(agent, {'fixed_uv', 'fixed_uv_discount'})
    if useastruct
        p.alpha  = params(strcmpi(astruct.parnames, 'alpha'));
        p.tau  = params(strcmpi(astruct.parnames, 'tau'));
    else
        p.beta = params(2);
        p.alpha = params(3);
        p.tau = params(4);
    end
elseif strcmpi(agent, 'null')
   %no parameters
end

if isfield(p, 'prop_spread') && length(p.prop_spread) > 0 && (p.prop_spread < 0 || p.prop_spread > 1), error('prop_spread outside of bounds'); end

end