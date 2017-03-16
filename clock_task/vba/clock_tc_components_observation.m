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
%       (last) sticky (sticky choice influence that evolves as a function of weighted RT history) [optional]
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

if ~isempty(u), RT = u(1); end 
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
    %RTlocavg = hidden(9)*inG.RTrescale; %not used in observation function
    expsign = hidden(10);
end

gx = 0; %empty predicted RT
phipos = 1; %current position within phi vector

%handle K + lambda
if ~isempty(regexp(inG.tcvariant, '^K_?', 'once')) %for K only model, there is no underscore after K
    K = unifinv(fastnormcdf(phi(phipos)), 0, inG.maxRT);
    gx = gx + K;
    phipos = phipos + 1;
end

if ~isempty(regexp(inG.tcvariant, '_Sticky', 'once'))
  %sticky choice is a function of two 0..1 distributed parameters
  %lambda and decay: lambda * sticky(t). sticky(t) = RT(t-1) + d*sticky(t-1)
  lambda = 1 ./ (1+exp(-phi(phipos))); %exponential transform to 0..1

  %sticky should be the last hidden state
  sticky = hidden(length(hidden));
  gx = gx + lambda * sticky;
  phipos = phipos + 1;
end

if ~isempty(regexp(inG.tcvariant, '_Lambda', 'once'))
    lambda = 1 ./ (1+exp(-phi(phipos))); %exponential transform to 0..1
    gx = gx + lambda*RT_prev;
    phipos = phipos + 1;
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
    nu = inG.maxNu ./ (1+exp(-phi(phipos))); %exponential transform to 0..10
    gx = gx + nu*(bestRT - meanRT);
    phipos = phipos + 1;
end

if ~isempty(regexp(inG.tcvariant, '_Rho|_Epsilon', 'once'))
    
    %compute the means and variances of fast and slow betas
    mean_slow = a_slow/(a_slow + b_slow);
    mean_fast = a_fast/(a_fast + b_fast);
    var_slow = a_slow * b_slow / ( (a_slow + b_slow)^2 * (a_slow + b_slow + 1) );
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
        rho = inG.rhoMultiply * gaminv(fastnormcdf(phi(phipos)), 2, 2); %use precompiled std norm cdf code for speed

        phipos = phipos + 1;
        gx = gx + rho*(mean_slow - mean_fast);
    end
    
    if ~isempty(regexp(inG.tcvariant, '_Epsilon', 'once'))
        %need a long-tailed gamma approximation? or just allow it to be negative and estimate accordingly?
        %for now, keep a non-negative number (otherwise need to implement sticky choice)
        
        %epsilon = inG.epsilonMultiply * gaminv(cdf('Normal', phi(phipos), inG.priors.muPhi(phipos), inG.priors.SigmaPhi(phipos,phipos)), 2, 2);
        
        %uniform 0..X variant
        %version allowing for mean and variance
        %epsilon = unifinv(normcdf(phi(phipos), inG.priors.muPhi(phipos), inG.priors.SigmaPhi(phipos,phipos)), 0, inG.maxEpsilon);
        
        %standard normal -> uniform variant (faster since it uses compiled code)
        %epsilon = unifinv(fastnormcdf(phi(phipos)), 0, inG.maxEpsilon);
        
        %try exponential distribution approach with mean of 1000 to yield a 0..~10000 distribution
        epsilon = expinv(fastnormcdf(phi(phipos)), inG.expEpsilonMean);

        %gaussian epsilon, allowing for negative (uncertainty averse)
        %epsilon = phi(phipos)*inG.epsilonMultiply;
        phipos = phipos + 1;

        gx = gx + epsilon*explore;
    end
end

if ~isempty(regexp(inG.tcvariant, '_Multinomial', 'once'))
  %gx should already contain the predicted RT (expectation)
  %now we just need to add a spread around it
  sd_U = (sqrt(var_slow) + sqrt(var_fast)) / 2; %use uncertainty to scale spread
  
  %compute SD of Gaussian
  %spread = phi(phipos)*sd_U*10; %x10 scale-up
  %spread = phi(phipos)*10; %x10 scale-up
  spread = phi(phipos);
  %spread = 4;
  
  tvec = 1:inG.ntimesteps; %vector of timesteps at which to evaluate gaussian

  %not using an untruncated basis function for normalization since that was a feature to scale height to 1.0, whereas here we want AUC=1.0 (PDF)
  %refspread = sum(gaussmf(min(tvec)-range(tvec):max(tvec)+range(tvec), [spread, median(tvec)])); %AUC of Gaussian in middle of interval

  %put predicted RT onto timestep grid
  gxrnd = round(gx*inG.ntimesteps/inG.range_RT);
  if gxrnd < 0, gxrnd=1; end %handle 0 timestep problem
  if gxrnd > inG.ntimesteps, gxrnd=inG.ntimesteps; end %handle overprediction problem
  
  %compute gaussian spread function around gxrnd (predicted RT timebin)
  elig = gaussmf(tvec, [spread, gxrnd]);

  %compute sum of area under the curve of the gaussian function
  auc=sum(elig);

  %divide gaussian update function by its sum so that AUC=1.0 (PDF)
  p_choice = elig/auc;
  
  %no softmax at the moment
  %p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature
  
  %return the normalized Gaussian as probabilities of each timestep (since they sum to 1)
  gx = p_choice;
end
