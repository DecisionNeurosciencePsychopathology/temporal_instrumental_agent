function  [fx] = clock_tc_components_evolution(hidden, theta, u, inF)
% evolution function for Frank TC model
%
% IN:
%   - hidden ("hidden states" to be updated) :
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
%   - theta (evolution parameters) : 
%       (1) alphaG (lr for positive PEs)
%       (2) alphaN (lr for negative PEs)
%
%   - u (inputs) :
%       (1) RT (RT on trial t)
%       (2) RT_prev (RT[t-1] -- used in observation function for resetting explore)
%       (3) reward (obtained on trial t)
%       (4) rew_max (maximum reward up to trial t -- used for tracking bestRT)
%       (5) rew_std (reward standard deviation up to trial t -- used for tracking bestRT)
%
%   - inF : struct of input options
%
% OUT:
%   - fx: evolved hidden states (in same order as input)

%initialize dummy variables (overridden if model uses them)
V = 0; reward = 0; rew_max = 0; rew_std = 0; bestRT = 0; Go = 0; NoGo = 0; alphaG = 0; alphaN = 0;

if ~isempty(regexp(inF.tcvariant, '_Nu', 'once'))
    bestRT = hidden(1)*inF.RTrescale;
    V = hidden(2);
end

if ~isempty(regexp(inF.tcvariant, '_AlphaG', 'once'))
    Go = hidden(3);
    alphaG = inF.maxAlpha ./ (1+exp(-theta(1))); %use exponential transform to rescale 0..1 and then multiply by 5 to be consistent with Frank (LR ranges from 0..5)
end

if ~isempty(regexp(inF.tcvariant, '_AlphaN', 'once'))
    NoGo = hidden(4);
    alphaN = inF.maxAlpha ./ (1+exp(-theta(2)));
end

if ~isempty(regexp(inF.tcvariant, '_Rho|_Epsilon', 'once'))
    a_fast = hidden(5);
    b_fast = hidden(6);
    a_slow = hidden(7);
    b_slow = hidden(8);
    RTlocavg = hidden(9)*inF.RTrescale;
    %expsign = hidden(10); %value not used in evolution (just set to +1/-1)
end

if length(u) > 0, RT = u(1); end 
%RT_prev is u(2); used in observation function
if length(u) > 2, reward = u(3); end
if length(u) > 3, rew_max = u(4); end
if length(u) > 4, rew_std = u(5); end

alphaV = 0.1; %Fixed as in Frank
V_new = V + alphaV*(reward - V); %critic value

bestRT_new = bestRT; %leave unchanged unless criterion below is met
if reward > V && reward >= (rew_max - rew_std)
    bestRT_new = RT;
end

%Go and NoGo terms only depend on reward and expected value.. seems like there should be a way to avoid their being hidden states..
Go_new = Go; %carry forward prior estimate unless updated by PPE
NoGo_new = NoGo;
if reward > V
    Go_new = Go - alphaG*(V - reward);
elseif reward <= V
    NoGo_new = NoGo + alphaN*(V - reward);
end

%process updates to beta distributions (for rho and epsilon models)
if ~isempty(regexp(inF.tcvariant, '_Epsilon|_Rho', 'once'))    
    a_slow_new = a_slow; b_slow_new = b_slow; %carry forward beta parameters unless updated below
    a_fast_new = a_fast; b_fast_new = b_fast;
    if RT > RTlocavg %update slow beta
        expsign_new = -1; %will speed RT in observation function if more uncertain about fast responses
        if reward > V %PE+
            a_slow_new = a_slow + 1;
        else %PE-
            b_slow_new = b_slow + 1;
        end
    elseif RT < RTlocavg && RT > 0 %update fast beta (verify that RT != 0 since that indicates non-response)
        expsign_new = +1; %will slow RT in observation function if more uncertain about slow responses
        if reward > V %PE+
            a_fast_new = a_fast + 1;
        else %PE-
            b_fast_new = b_fast + 1;
        end
    end
    
    %update estimate of recent/"learned" RTlocavg
    RTlocavg_new = RTlocavg + alphaV*(reward - RTlocavg);
end

fx = []; 
if ~isempty(regexp(inF.tcvariant, '_Nu', 'once'))
    fx = [fx; bestRT_new/inF.RTrescale; V_new];
end

if ~isempty(regexp(inF.tcvariant, '_AlphaG', 'once'))
    fx = [fx; Go_new];
end

if ~isempty(regexp(inF.tcvariant, '_AlphaN', 'once'))
    fx = [fx; NoGo_new];
end

if ~isempty(regexp(inF.tcvariant, '_Rho|_Epsilon'))
    fx = [fx; a_fast_new; b_fast_new; a_slow_new; b_slow_new; RTlocavg_new/inF.RTrescale; expsign_new];
end