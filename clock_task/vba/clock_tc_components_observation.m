function  [ gx ] = clock_tc_components_observation(hidden,phi,u,inG)
% Frank TC observation (choice) function

% INPUT
%   - hidden : evolved hidden states from evolution function
%       (1)  bestRT (RT associated with largest reward up to this trial)
%       (2)  V (expected value of choice)
%       (3)  Go (speeding of RT by PPE)
%       (4)  NoGo (slowing of RT by NPE)
%       (5)  a_fast (alpha parameter governing beta distribution for fast RTs)
%       (6)  b_fast (beta parameter governing beta distribution for fast RTs)
%       (7)  a_slow  (alpha parameter governing beta distribution for slow RTs)
%       (8)  b_slow (beta parameter governing beta distribution for slow RTs)
%       (9)  RTlocavg (local average/"learned" RT)
%       (10) expsign (+1/-1 scalar indicating whether uncertainty-driven exploration speeds or slows RT
%
%   - phi : observation parameters
%       (1) K (intercept)
%       (2) lambda (autocorrelation)
%       (3) nu (go for the gold)
%       (4) rho (mean value difference for fast/slow)
%       (5) epsilon (explore based on rel. diff in fast/slow uncertainty)
%       
%   - u (inputs) :
%       (1) RT (RT on trial t)
%       (2) RT_prev (RT[t-1] -- used in observation function for resetting explore)
%       (3) reward (obtained on trial t)
%       (4) rew_max (maximum reward up to trial t -- used for tracking bestRT)
%       (5) rew_std (reward standard deviation up to trial t -- used for tracking bestRT)
%
% - inG : input structure
%
% OUTPUT
% - gx : p(chosen|hidden)

%pull out input parameters

if length(u) > 0, RT = u(1); end 
if length(u) > 1, RT_prev = u(2); end
    
%these are not used in observation function
%if length(u) > 2, reward = u(3); end
%if length(u) > 3, rew_max = u(4); end
%if length(u) > 4, rew_std = u(5); end


%pull out hidden states
if ~isempty(regexp(inG.tcvariant, '_Nu', 'once'))
    bestRT = hidden(1)*inG.RTrescale;
    %V = hidden(2); %not used in observation
end

if ~isempty(regexp(inG.tcvariant, '_AlphaG', 'once'))
    Go = hidden(3);
end

if ~isempty(regexp(inG.tcvariant, '_AlphaN', 'once'))
    NoGo = hidden(4);
end

if ~isempty(regexp(inG.tcvariant, '_Rho|_Epsilon', 'once'))
    a_fast = hidden(5);
    b_fast = hidden(6);
    a_slow = hidden(7);
    b_slow = hidden(8);
    %RTlocavg = hidden(9)*inG.RTrescale;
    expsign = hidden(10);
end

gx = 0; %empty predicted RT

%handle K + lambda
if ~isempty(regexp(inG.tcvariant, 'K_', 'once'))
    K = unifinv(fastnormcdf(phi(1)), 0, inG.maxRT);
    %K = unifinv(interpn(inG.stdnormgrid, inG.stdnormcdf, phi(1), 'linear'), 0, inG.maxRT);
    gx = gx + K;
end

if ~isempty(regexp(inG.tcvariant, '_Lambda', 'once'))
    lambda = 1 ./ (1+exp(-phi(2))); %exponential transform to 0..1
    gx = gx + lambda*RT_prev;
end

if ~isempty(regexp(inG.tcvariant, '_AlphaG', 'once'))
    gx = gx - Go;
end

if ~isempty(regexp(inG.tcvariant, '_AlphaN', 'once'))
    gx = gx + NoGo;
end

if ~isempty(regexp(inG.tcvariant, '_Nu', 'once'))
    %need overall average RT for the block for "Go for the gold"
    meanRT = inG.meanRT;
    nu = inG.maxNu ./ (1+exp(-phi(3))); %exponential transform to 0..100
    gx = gx + nu*(bestRT - meanRT);
end

if ~isempty(regexp(inG.tcvariant, '_Rho|_Epsilon', 'once'))
    
    %compute the means and variances of fast and slow betas
    mean_slow = a_slow/(a_slow + b_slow);
    mean_fast = a_fast/(a_fast + b_fast);
    var_slow = a_slow*b_slow/( (a_slow + b_slow)^2 * (a_slow + b_slow + 1) );
    var_fast = a_fast * b_fast / ( (a_fast + b_fast)^2 * (a_fast + b_fast + 1) );
    
    explore=0; %default to 0 explore term unless updated in tc_evolution
    if expsign == 1 %RT was faster than avg, so potentially slow down
        explore = expsign*(sqrt(var_slow) - sqrt(var_fast));
    elseif expsign == -1 %RT was slower than average, so potentially speed up
        explore = expsign*(sqrt(var_fast) - sqrt(var_slow));
    end
    
    %nullify exploration parameter if the shift in RT from t-1 to t was already
    %in the same direction as the explore parameter would achieve
    if RT < RT_prev && explore < 0
        explore=0;
    elseif RT > RT_prev && explore > 0
        explore=0;
    end
    
    if ~isempty(regexp(inG.tcvariant, '_Rho', 'once'))        
        %use a long-tailed gamma(2,2) distribution to approximate 0..10000 parameter
        %at probability .999, gamma is ~18, so need to multiply by about 400 to achieve 0..10000 scaling
        %rho = inG.rhoMultiply * gaminv(normcdf(phi(4), inG.priors.muPhi(4), inG.priors.SigmaPhi(4,4)), 2, 2);
        rho = inG.rhoMultiply * gaminv(fastnormcdf(phi(4)), 2, 2); %use precompiled std norm cdf code for speed
        %rho = inG.rhoMultiply * gaminv(interpn(inG.stdnormgrid, inG.stdnormcdf, phi(4), 'linear'), 2, 2);

        gx = gx + rho*(mean_slow - mean_fast);
    end
    
    if ~isempty(regexp(inG.tcvariant, '_Epsilon', 'once'))
        %need a long-tailed gamma approximation? or just allow it to be negative and estimate accordingly?
        %for now, keep a non-negative number (otherwise need to implement sticky choice)
        
        %epsilon = inG.epsilonMultiply * gaminv(cdf('Normal', phi(5), inG.priors.muPhi(5), inG.priors.SigmaPhi(5,5)), 2, 2);
        %uniform 0..80000 variant
        %epsilon = unifinv(normcdf(phi(5), inG.priors.muPhi(5), inG.priors.SigmaPhi(5,5)), 0, 80000);
        epsilon = unifinv(fastnormcdf(phi(5)), 0, 80000);
        %epsilon = unifinv(interpn(inG.stdnormgrid, inG.stdnormcdf, phi(5), 'linear'), 0, 80000);

        gx = gx + epsilon*explore;
    end
end