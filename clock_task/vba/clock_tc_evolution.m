function  [fx] = clock_tc_evolution(hidden, theta, u, inF)
% evolution function for Frank TC model
%
% IN:
%   - hidden ("hidden states" to be updated) :
%       (1)  bestRT (RT associated with largest reward up to this trial)
%       (2)  a_fast (alpha parameter governing beta distribution for fast RTs)
%       (3)  b_fast (beta parameter governing beta distribution for fast RTs)
%       (4)  a_slow  (alpha parameter governing beta distribution for slow RTs)
%       (5)  b_slow (beta parameter governing beta distribution for slow RTs)
%       (6)  V (expected value of choice)
%       (7)  Go (speeding of RT by PPE)
%       (8)  NoGo (slowing of RT by NPE)
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

%extract elements of input
RT = u(1); %current RT
%RT_prev = u(2); %previous RT (not used here in evolution function)
reward = u(3); %reward obtained on trial t
rew_max = u(4); %max reward up to trial t
rew_std = u(5); %reward standard deviation up to trial t

alphaG = inF.maxAlpha ./ (1+exp(-theta(1))); %use exponential transform to rescale 0..1 and then multiply by 5 to be consistent with Frank (LR ranges from 0..5)
alphaN = inF.maxAlpha ./ (1+exp(-theta(2)));

%any model component that changes with trial is treated as a hidden state.
%any model component that doesn't evolve is placed in the evolution function.

%extract hidden states to update
bestRT = hidden(1)*inF.RTrescale;
a_fast = hidden(2);
b_fast = hidden(3);
a_slow = hidden(4);
b_slow = hidden(5);
V = hidden(6);
Go = hidden(7);
NoGo = hidden(8);
RTlocavg = hidden(9)*inF.RTrescale;
%expsign = hidden(10);

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

%process updates to beta distributions
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
RTlocavg_new = RTlocavg + alphaV*(RT - RTlocavg);

fx = [bestRT_new/inF.RTrescale; a_fast_new; b_fast_new; a_slow_new; b_slow_new; ...
    V_new; Go_new; NoGo_new; RTlocavg_new/inF.RTrescale; expsign_new];
