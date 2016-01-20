function  [fx] = clock_tc_evolution(hidden, theta, u, inF)
% evolution function for K + lambda model
%
% IN:
%   - hidden ("hidden states" to be updated) : NONE
%
%   - theta (evolution parameters) : NONE
%
%   - u (inputs) :
%       (1) RT (RT on trial t)
%       (2) RT_prev (RT[t-1] -- used in observation function for resetting explore)
%       (3) reward (obtained on trial t)
%
%   - inF : struct of input options
%
% OUT:
%   - fx: evolved hidden states (in same order as input)

%extract elements of input
%RT = u(1); %current RT
%RT_prev = u(2); %previous RT (not used here in evolution function)

%empty evolution function because there are no hidden states

fx = [];